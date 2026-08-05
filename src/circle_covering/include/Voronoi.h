#pragma once
// STL
#include <limits>
#include <set>
#include <memory>
#include <unordered_set>
#include <unordered_map>
// My includes
#include "Clip.h"

// Box
class Box
{
public:
    // Be careful, y-axis is oriented to the top like in math
    enum class Side : int
    {
        LEFT,
        BOTTOM,
        RIGHT,
        TOP
    };

    struct Intersection
    {
        Side side;
        Point point;
    };

    double left;
    double bottom;
    double right;
    double top;

    bool contains(const Point &point) const;
    Intersection getFirstIntersection(const Point &origin, const Point &direction) const;                                  // Useful for Fortune's algorithm
    int getIntersections(const Point &origin, const Point &destination, std::array<Intersection, 2> &intersections) const; // Useful for diagram intersection

private:
    static constexpr double EPSILON = 1e-7;
};

// VoronoiDiagram
class FortuneAlgorithm;
class VoronoiDiagram
{
public:
    struct HalfEdge;
    struct Cell;

    struct Site
    {
        std::size_t index;
        Point point;
        Cell *cell;
    };

    struct Vertex
    {
        Point point;

    private:
        friend VoronoiDiagram;
        std::list<Vertex>::iterator it;
    };

    struct HalfEdge
    {
        Vertex *origin = nullptr;
        Vertex *destination = nullptr;
        HalfEdge *twin = nullptr;
        Cell *incidentFace;
        HalfEdge *prev = nullptr;
        HalfEdge *next = nullptr;

    private:
        friend VoronoiDiagram;
        std::list<HalfEdge>::iterator it;
    };

    struct Cell
    {
        Site *site;
        HalfEdge *outerComponent;
    };

    VoronoiDiagram(const std::vector<Point> &points);

    // Remove copy operations
    VoronoiDiagram(const VoronoiDiagram &) = delete;
    VoronoiDiagram &operator=(const VoronoiDiagram &) = delete;

    // Move operations
    VoronoiDiagram(VoronoiDiagram &&) = default;
    VoronoiDiagram &operator=(VoronoiDiagram &&) = default;

    // Accessors
    const std::vector<Site> &getSites() const;
    Site *getSite(std::size_t i);
    std::size_t getNbSites() const;
    Cell *getCell(std::size_t i);
    const std::list<Vertex> &getVertices() const;
    const std::vector<Point> getPoints(std::size_t i) const;
    const std::list<HalfEdge> &getHalfEdges() const;

    // Intersection with a box
    bool intersect(Box box);
    std::vector<std::vector<Point>> result();
    std::vector<std::vector<Point>> result_with_neighbors();

    std::vector<Point> result(size_t);

private:
    std::vector<Site> mSites;
    std::vector<Cell> mCells;
    std::list<Vertex> mVertices;
    std::list<HalfEdge> mHalfEdges;

    // Diagram construction
    friend FortuneAlgorithm;

    Vertex *createVertex(Point point);
    Vertex *createCorner(Box box, Box::Side side);
    HalfEdge *createHalfEdge(Cell *cell);

    // Intersection with a box
    void link(Box box, HalfEdge *start, Box::Side startSide, HalfEdge *end, Box::Side endSide);
    void removeVertex(Vertex *vertex);
    void removeHalfEdge(HalfEdge *halfEdge);
};

// Arc
class Event;
struct Arc
{
    enum class Color
    {
        RED,
        BLACK
    };
    enum class Side
    {
        Left,
        Right
    };

    // Hierarchy
    Arc *parent;
    Arc *left;
    Arc *right;
    // Diagram
    VoronoiDiagram::Site *site;
    VoronoiDiagram::HalfEdge *leftHalfEdge;
    VoronoiDiagram::HalfEdge *rightHalfEdge;
    Event *event;
    // Optimizations
    Arc *prev;
    Arc *next;
    // Only for balancing
    Color color;
    // To know if the arc is towards -inf or +inf
    Side side;
};

// Event
class Event
{
public:
    enum class Type
    {
        SITE,
        CIRCLE
    };

    // Site event
    Event(VoronoiDiagram::Site *site);
    // Circle event
    Event(double y, Point point, Arc *arc);

    const Type type;
    double x;
    double y;
    int index;
    // Site event
    VoronoiDiagram::Site *site;
    // Circle event
    Point point;
    Arc *arc;
};

bool operator<(const Event &lhs, const Event &rhs);

// Beachline
class Beachline
{
public:
    Beachline();
    ~Beachline();

    // Remove copy operations

    Beachline(const Beachline &) = delete;
    Beachline &operator=(const Beachline &) = delete;

    // Move operations
    Beachline(Beachline &&other);
    Beachline &operator=(Beachline &&other);
    Arc *createArc(VoronoiDiagram::Site *site, Arc::Side side);
    bool isEmpty() const;
    bool isNil(const Arc *x) const;
    void setRoot(Arc *x);
    Arc *getLeftmostArc() const;
    Arc *locateArcAbove(const Point &point, double l) const;
    void insertBefore(Arc *x, Arc *y);
    void insertAfter(Arc *x, Arc *y);
    void replace(Arc *x, Arc *y);
    void remove(Arc *z);

private:
    Arc *mNil;
    Arc *mRoot;

    // Utility methods
    Arc *minimum(Arc *x) const;
    void transplant(Arc *u, Arc *v);

    // Fixup functions
    void insertFixup(Arc *z);
    void removeFixup(Arc *x);

    // Rotations
    void leftRotate(Arc *x);
    void rightRotate(Arc *y);
    double computeBreakpoint(const Point &point1, const Point &point2, double l, Arc::Side side) const;
    void free(Arc *x);
};

// PriorityQueue
template <typename T>
class PriorityQueue
{
public:
    PriorityQueue()
    {
    }

    // Accessors

    bool isEmpty() const
    {
        return mElements.empty();
    }

    // Operations

    std::unique_ptr<T> pop()
    {
        swap(0, mElements.size() - 1);
        auto top = std::move(mElements.back());
        mElements.pop_back();
        siftDown(0);
        return top;
    }

    void push(std::unique_ptr<T> elem)
    {
        elem->index = mElements.size();
        mElements.emplace_back(std::move(elem));
        siftUp(mElements.size() - 1);
    }

    void update(std::size_t i)
    {
        int parent = getParent(i);
        if (parent >= 0 && *mElements[parent] < *mElements[i])
            siftUp(i);
        else
            siftDown(i);
    }

    void remove(std::size_t i)
    {
        swap(i, mElements.size() - 1);
        mElements.pop_back();
        if (i < mElements.size())
            update(i);
    }

    // Print

    std::ostream &print(std::ostream &os, int i = 0, std::string tabs = "") const
    {
        if (i < mElements.size())
        {
            os << tabs << *mElements[i] << std::endl;
            display(getLeftChild(i), tabs + '\t');
            display(getRightChild(i), tabs + '\t');
        }
        return os;
    }

private:
    std::vector<std::unique_ptr<T>> mElements;

    // Accessors

    int getParent(int i) const
    {
        return (i + 1) / 2 - 1;
    }

    std::size_t getLeftChild(std::size_t i) const
    {
        return 2 * (i + 1) - 1;
    }

    std::size_t getRightChild(std::size_t i) const
    {
        return 2 * (i + 1);
    }

    // Operations

    void siftDown(std::size_t i)
    {
        std::size_t left = getLeftChild(i);
        std::size_t right = getRightChild(i);
        std::size_t j = i;
        if (left < mElements.size() && *mElements[j] < *mElements[left])
            j = left;
        if (right < mElements.size() && *mElements[j] < *mElements[right])
            j = right;
        if (j != i)
        {
            swap(i, j);
            siftDown(j);
        }
    }

    void siftUp(std::size_t i)
    {
        int parent = getParent(i);
        if (parent >= 0 && *mElements[parent] < *mElements[i])
        {
            swap(i, parent);
            siftUp(parent);
        }
    }

    inline void swap(std::size_t i, std::size_t j)
    {
        std::swap(mElements[i], mElements[j]);
        mElements[i]->index = i;
        mElements[j]->index = j;
    }
};

// FortuneAlgorithm
class FortuneAlgorithm
{
public:
    FortuneAlgorithm(std::vector<Point> points);
    ~FortuneAlgorithm();

    void construct();
    bool bound(Box box);

    VoronoiDiagram getDiagram();

private:
    VoronoiDiagram mDiagram;
    Beachline mBeachline;
    PriorityQueue<Event> mEvents;
    double mBeachlineY;

    // Algorithm
    void handleSiteEvent(Event *event);
    void handleCircleEvent(Event *event);

    // Arcs
    Arc *breakArc(Arc *arc, VoronoiDiagram::Site *site);
    void removeArc(Arc *arc, VoronoiDiagram::Vertex *vertex);

    // Breakpoint
    bool isMovingRight(const Arc *left, const Arc *right) const;
    double getInitialX(const Arc *left, const Arc *right, bool movingRight) const;

    // Edges
    void addEdge(Arc *left, Arc *right);
    void setOrigin(Arc *left, Arc *right, VoronoiDiagram::Vertex *vertex);
    void setDestination(Arc *left, Arc *right, VoronoiDiagram::Vertex *vertex);
    void setPrevHalfEdge(VoronoiDiagram::HalfEdge *prev, VoronoiDiagram::HalfEdge *next);

    // Events
    void addEvent(Arc *left, Arc *middle, Arc *right);
    void deleteEvent(Arc *arc);
    Point computeConvergencePoint(const Point &point1, const Point &point2, const Point &point3, double &y) const;

    // Bounding

    struct LinkedVertex
    {
        VoronoiDiagram::HalfEdge *prevHalfEdge;
        VoronoiDiagram::Vertex *vertex;
        VoronoiDiagram::HalfEdge *nextHalfEdge;
    };
};
