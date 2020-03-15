#ifndef CAVC_STATICSPATIALINDEX_H
#define CAVC_STATICSPATIALINDEX_H
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <stack>
#include <vector>

namespace cavc {
template <typename Real, std::size_t NodeSize = 16> class StaticSpatialIndex {
public:
  StaticSpatialIndex(std::size_t numItems) {
    assert(numItems > 0 && "number of items must be greater than 0");
    static_assert(NodeSize >= 2 && NodeSize <= 65535, "node size must be between 2 and 65535");
    // calculate the total number of nodes in the R-tree to allocate space for
    // and the index of each tree level (used in search later)
    m_numItems = numItems;
    std::size_t n = numItems;
    std::size_t numNodes = numItems;
    m_levelBounds.push_back(n * 4);
    do {
      n = static_cast<std::size_t>(std::ceil(static_cast<float>(n) / NodeSize));
      numNodes += n;
      m_levelBounds.push_back(numNodes * 4);
    } while (n != 1);

    m_boxes.resize(numNodes * 4);
    m_indices.resize(numNodes);
    m_pos = 0;
    m_minX = std::numeric_limits<Real>::infinity();
    m_minY = std::numeric_limits<Real>::infinity();
    m_maxX = -1 * std::numeric_limits<Real>::infinity();
    m_maxY = -1 * std::numeric_limits<Real>::infinity();
  }

  Real MinX() const { return m_minX; }
  Real MinY() const { return m_minY; }
  Real MaxX() const { return m_maxX; }
  Real MaxY() const { return m_maxY; }

  void add(Real minX, Real minY, Real maxX, Real maxY) {
    std::size_t index = m_pos >> 2;
    m_indices[index] = index;
    m_boxes[m_pos++] = minX;
    m_boxes[m_pos++] = minY;
    m_boxes[m_pos++] = maxX;
    m_boxes[m_pos++] = maxY;

    if (minX < m_minX)
      m_minX = minX;
    if (minY < m_minY)
      m_minY = minY;
    if (maxX > m_maxX)
      m_maxX = maxX;
    if (maxY > m_maxY)
      m_maxY = maxY;
  }

  void finish() {
    assert(m_pos >> 2 == m_numItems && "added item count should equal static size given");
    Real width = m_maxX - m_minX;
    Real height = m_maxY - m_minY;
    std::vector<std::uint32_t> hilbertValues(m_numItems);

    std::size_t pos = 0;

    for (std::size_t i = 0; i < m_numItems; ++i) {
      pos = 4 * i;
      Real minX = m_boxes[pos++];
      Real minY = m_boxes[pos++];
      Real maxX = m_boxes[pos++];
      Real maxY = m_boxes[pos++];

      // hilbert max input value for x and y
      const Real hilbertMax = static_cast<Real>((1 << 16) - 1);
      // mapping the x and y coordinates of the center of the box to values in the range
      // [0 -> n - 1] such that the min of the entire set of bounding boxes maps to 0 and the max of
      // the entire set of bounding boxes maps to n - 1 our 2d space is x: [0 -> n-1] and
      // y: [0 -> n-1], our 1d hilbert curve value space is d: [0 -> n^2 - 1]
      Real x = std::floor(hilbertMax * ((minX + maxX) / 2 - m_minX) / width);
      std::uint32_t hx = static_cast<std::uint32_t>(x);
      Real y = std::floor(hilbertMax * ((minY + maxY) / 2 - m_minY) / height);
      std::uint32_t hy = static_cast<std::uint32_t>(y);
      hilbertValues[i] = hilbertXYToIndex(hx, hy);
    }

    // sort items by their Hilbert value (for packing later)
    sort(hilbertValues, m_boxes, m_indices, 0, m_numItems - 1);

    // generate nodes at each tree level, bottom-up
    pos = 0;
    for (std::size_t i = 0; i < m_levelBounds.size() - 1; i++) {
      auto end = m_levelBounds[i];

      // generate a parent node for each block of consecutive <nodeSize> nodes
      while (pos < end) {
        auto nodeMinX = std::numeric_limits<Real>::infinity();
        auto nodeMinY = std::numeric_limits<Real>::infinity();
        auto nodeMaxX = -1 * std::numeric_limits<Real>::infinity();
        auto nodeMaxY = -1 * std::numeric_limits<Real>::infinity();
        auto nodeIndex = pos;

        // calculate bbox for the new node
        for (std::size_t j = 0; j < NodeSize && pos < end; j++) {
          auto minX = m_boxes[pos++];
          auto minY = m_boxes[pos++];
          auto maxX = m_boxes[pos++];
          auto maxY = m_boxes[pos++];
          if (minX < nodeMinX)
            nodeMinX = minX;
          if (minY < nodeMinY)
            nodeMinY = minY;
          if (maxX > nodeMaxX)
            nodeMaxX = maxX;
          if (maxY > nodeMaxY)
            nodeMaxY = maxY;
        }

        // add the new node to the tree data
        m_indices[m_pos >> 2] = nodeIndex;
        m_boxes[m_pos++] = nodeMinX;
        m_boxes[m_pos++] = nodeMinY;
        m_boxes[m_pos++] = nodeMaxX;
        m_boxes[m_pos++] = nodeMaxY;
      }
    }
  }

  // Visit all the bounding boxes in the spatial index. Visitor function has the signature
  // void(Real xmin, Real ymin, Real xmax, Real ymax, std::size_t level).
  template <typename F> void visitBoundingBoxes(F &&visitor) const {
    std::size_t nodeIndex = m_boxes.size() - 4;
    std::size_t level = m_levelBounds.size() - 1;

    std::vector<std::size_t> stack;
    stack.reserve(16);

    bool done = false;
    while (!done) {
      auto end = std::min(nodeIndex + NodeSize * 4, m_levelBounds[level]);
      for (std::size_t pos = nodeIndex; pos < end; pos += 4) {
        auto index = m_indices[pos >> 2];
        visitor(m_boxes[pos], m_boxes[pos + 1], m_boxes[pos + 2], m_boxes[pos + 3], level);

        if (nodeIndex >= m_numItems * 4) {
          stack.push_back(index);
          stack.push_back(level - 1);
        }
      }

      if (stack.size() > 1) {
        level = stack.back();
        stack.pop_back();
        nodeIndex = stack.back();
        stack.pop_back();
      } else {
        done = true;
      }
    }
  }

  // Query the spatial index adding indexes to the results vector given.
  void query(Real minX, Real minY, Real maxX, Real maxY, std::vector<std::size_t> &results) const {
    auto visitor = [&](std::size_t index) {
      results.push_back(index);
      return true;
    };

    visitQuery(minX, minY, maxX, maxY, visitor);
  }

  // Query the spatial index, invoking a visitor function for each index that overlaps the bounding
  // box given. Visitor function has the signature bool(std::size_t index), if visitor returns false
  // the query stops early, otherwise the query continues.
  template <typename F>
  void visitQuery(Real minX, Real minY, Real maxX, Real maxY, F &&visitor) const {
    assert(m_pos == m_boxes.size() && "data not yet indexed - call Finish() before querying");

    auto nodeIndex = m_boxes.size() - 4;
    auto level = m_levelBounds.size() - 1;

    // stack for traversing nodes
    std::vector<std::size_t> stack;
    // reserve some space to avoid repeated small allocations
    stack.reserve(16);

    auto done = false;

    while (!done) {
      // find the end index of the node
      auto end = std::min(nodeIndex + NodeSize * 4, m_levelBounds[level]);

      // search through child nodes
      for (std::size_t pos = nodeIndex; pos < end; pos += 4) {
        auto index = m_indices[pos >> 2];
        // check if node bbox intersects with query bbox
        if (maxX < m_boxes[pos])
          continue; // maxX < nodeMinX
        if (maxY < m_boxes[pos + 1])
          continue; // maxY < nodeMinY
        if (minX > m_boxes[pos + 2])
          continue; // minX > nodeMaxX
        if (minY > m_boxes[pos + 3])
          continue; // minY > nodeMaxY

        if (nodeIndex < m_numItems * 4) {
          done = !visitor(index);
          if (done) {
            break;
          }
        } else {
          // push node index and level for further traversal
          stack.push_back(index);
          stack.push_back(level - 1);
        }
      }

      if (stack.size() > 1) {
        level = stack.back();
        stack.pop_back();
        nodeIndex = stack.back();
        stack.pop_back();
      } else {
        done = true;
      }
    }
  }

  static std::uint32_t hilbertXYToIndex(std::uint32_t x, std::uint32_t y) {
    std::uint32_t a = x ^ y;
    std::uint32_t b = 0xFFFF ^ a;
    std::uint32_t c = 0xFFFF ^ (x | y);
    std::uint32_t d = x & (y ^ 0xFFFF);

    std::uint32_t A = a | (b >> 1);
    std::uint32_t B = (a >> 1) ^ a;
    std::uint32_t C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
    std::uint32_t D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;

    a = A;
    b = B;
    c = C;
    d = D;
    A = ((a & (a >> 2)) ^ (b & (b >> 2)));
    B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));
    C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
    D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));

    a = A;
    b = B;
    c = C;
    d = D;
    A = ((a & (a >> 4)) ^ (b & (b >> 4)));
    B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));
    C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
    D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));

    a = A;
    b = B;
    c = C;
    d = D;
    C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
    D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));

    a = C ^ (C >> 1);
    b = D ^ (D >> 1);

    std::uint32_t i0 = x ^ y;
    std::uint32_t i1 = b | (0xFFFF ^ (i0 | a));

    i0 = (i0 | (i0 << 8)) & 0x00FF00FF;
    i0 = (i0 | (i0 << 4)) & 0x0F0F0F0F;
    i0 = (i0 | (i0 << 2)) & 0x33333333;
    i0 = (i0 | (i0 << 1)) & 0x55555555;

    i1 = (i1 | (i1 << 8)) & 0x00FF00FF;
    i1 = (i1 | (i1 << 4)) & 0x0F0F0F0F;
    i1 = (i1 | (i1 << 2)) & 0x33333333;
    i1 = (i1 | (i1 << 1)) & 0x55555555;

    return (i1 << 1) | i0;
  }

private:
  Real m_minX;
  Real m_minY;
  Real m_maxX;
  Real m_maxY;
  std::size_t m_numItems;
  std::vector<std::size_t> m_levelBounds;
  std::vector<Real> m_boxes;
  std::vector<std::size_t> m_indices;
  std::size_t m_pos;

  static void sort(std::vector<std::uint32_t> &values, std::vector<Real> &boxes,
                   std::vector<std::size_t> &indices, std::size_t left, std::size_t right) {

    if (left >= right)
      return;

    auto pivot = values[(left + right) >> 1];
    auto i = left - 1;
    auto j = right + 1;

    while (true) {
      do
        i++;
      while (values[i] < pivot);
      do
        j--;
      while (values[j] > pivot);
      if (i >= j)
        break;
      swap(values, boxes, indices, i, j);
    }

    sort(values, boxes, indices, left, j);
    sort(values, boxes, indices, j + 1, right);
  }

  static void swap(std::vector<std::uint32_t> &values, std::vector<Real> &boxes,
                   std::vector<std::size_t> &indices, std::size_t i, std::size_t j) {
    auto temp = values[i];
    values[i] = values[j];
    values[j] = temp;

    auto k = 4 * i;
    auto m = 4 * j;

    auto a = boxes[k];
    auto b = boxes[k + 1];
    auto c = boxes[k + 2];
    auto d = boxes[k + 3];
    boxes[k] = boxes[m];
    boxes[k + 1] = boxes[m + 1];
    boxes[k + 2] = boxes[m + 2];
    boxes[k + 3] = boxes[m + 3];
    boxes[m] = a;
    boxes[m + 1] = b;
    boxes[m + 2] = c;
    boxes[m + 3] = d;

    auto e = indices[i];
    indices[i] = indices[j];
    indices[j] = e;
  }
};
} // namespace cavc

#endif // CAVC_STATICSPATIALINDEX_H
