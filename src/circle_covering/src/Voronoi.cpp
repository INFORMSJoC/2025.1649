#include "Voronoi.h"

// GLOBAL
std::vector<std::vector<int>> neighbors;
extern double range;

// Box
bool Box::contains(const Point &point) const
{
    return (point.x >= left - EPSILON && point.x <= right + EPSILON) &&
           (point.y >= bottom - EPSILON && point.y <= top + EPSILON);
}

Box::Intersection Box::getFirstIntersection(const Point &origin, const Point &direction) const
{
    // origin must be in the box
    Intersection intersection;
    double t = std::numeric_limits<double>::infinity();
    if (direction.x > 0.0)
    {
        t = (right - origin.x) / direction.x;
        intersection.side = Side::RIGHT;
        intersection.point = origin + t * direction;
    }
    else if (direction.x < 0.0)
    {
        t = (left - origin.x) / direction.x;
        intersection.side = Side::LEFT;
        intersection.point = origin + t * direction;
    }
    if (direction.y > 0.0)
    {
        double newT = (top - origin.y) / direction.y;
        if (newT < t)
        {
            intersection.side = Side::TOP;
            intersection.point = origin + newT * direction;
        }
    }
    else if (direction.y < 0.0)
    {
        double newT = (bottom - origin.y) / direction.y;
        if (newT < t)
        {
            intersection.side = Side::BOTTOM;
            intersection.point = origin + newT * direction;
        }
    }
    return intersection;
}

int Box::getIntersections(const Point &origin, const Point &destination, std::array<Intersection, 2> &intersections) const
{
    // WARNING: If the intersection is a corner, both intersections are equals
    Point direction = destination - origin;
    std::array<double, 2> t;
    std::size_t i = 0; // index of the current intersection
    // Left
    if ((origin.x < left - EPSILON) || (destination.x < left - EPSILON))
    {
        t[i] = (left - origin.x) / direction.x;
        if ((t[i] > EPSILON) && (t[i] < 1.0 - EPSILON))
        {
            intersections[i].side = Side::LEFT;
            intersections[i].point = origin + t[i] * direction;
            if ((intersections[i].point.y >= bottom - EPSILON) && (intersections[i].point.y <= top + EPSILON))
                ++i;
        }
    }
    // Right
    if ((origin.x > right + EPSILON) || (destination.x > right + EPSILON))
    {
        t[i] = (right - origin.x) / direction.x;
        if ((t[i] > EPSILON) && (t[i] < 1.0 - EPSILON))
        {
            intersections[i].side = Side::RIGHT;
            intersections[i].point = origin + t[i] * direction;
            if ((intersections[i].point.y >= bottom - EPSILON) && (intersections[i].point.y <= top + EPSILON))
                ++i;
        }
    }
    // Bottom
    if ((origin.y < bottom - EPSILON) || (destination.y < bottom - EPSILON))
    {
        t[i] = (bottom - origin.y) / direction.y;
        if (i < 2 && ((t[i] > EPSILON) && (t[i] < 1.0 - EPSILON)))
        {
            intersections[i].side = Side::BOTTOM;
            intersections[i].point = origin + t[i] * direction;
            if ((intersections[i].point.x >= left - EPSILON) && (intersections[i].point.x <= right + EPSILON))
                ++i;
        }
    }
    // Top
    if ((origin.y > top + EPSILON) || (destination.y > top + EPSILON))
    {
        t[i] = (top - origin.y) / direction.y;
        if (i < 2 && ((t[i] > EPSILON) && (t[i] < 1.0 - EPSILON)))
        {
            intersections[i].side = Side::TOP;
            intersections[i].point = origin + t[i] * direction;
            if ((intersections[i].point.x >= left - EPSILON) && (intersections[i].point.x <= right + EPSILON))
                ++i;
        }
    }
    // Sort the intersections from the nearest to the farthest
    if (i == 2 && t[0] > t[1])
        std::swap(intersections[0], intersections[1]);
    return i;
}

// VoronoiDiagram
VoronoiDiagram::VoronoiDiagram(const std::vector<Point> &points)
{
    mSites.reserve(points.size());
    mCells.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); i++)
    {
        mSites.push_back(VoronoiDiagram::Site{i, points[i], nullptr});
        mCells.push_back(VoronoiDiagram::Cell{&mSites.back(), nullptr});
        mSites.back().cell = &mCells.back();
    }
}

const std::vector<VoronoiDiagram::Site> &VoronoiDiagram::getSites() const
{
    return mSites;
}

VoronoiDiagram::Site *VoronoiDiagram::getSite(std::size_t i)
{
    return &mSites[i];
}

std::size_t VoronoiDiagram::getNbSites() const
{
    return mSites.size();
}

VoronoiDiagram::Cell *VoronoiDiagram::getCell(std::size_t i)
{
    return &mCells[i];
}

const std::list<VoronoiDiagram::Vertex> &VoronoiDiagram::getVertices() const
{
    return mVertices;
}

const std::vector<Point> VoronoiDiagram::getPoints(std::size_t i) const
{
    std::vector<Point> vertices;
    Site site = mSites[i];
    Cell *cell = site.cell;
    HalfEdge *halfEdge = cell->outerComponent;
    HalfEdge *start = halfEdge;
    while (halfEdge != nullptr)
    {
        if (halfEdge->origin != nullptr && halfEdge->destination != nullptr)
        {
            Point origin = halfEdge->origin->point;
            vertices.push_back(origin);
        }
        halfEdge = halfEdge->next;
        if (halfEdge == start)
            break;
    }
    return vertices;
}

const std::list<VoronoiDiagram::HalfEdge> &VoronoiDiagram::getHalfEdges() const
{
    return mHalfEdges;
}

bool VoronoiDiagram::intersect(Box box)
{
    bool success = true;
    std::unordered_set<HalfEdge *> processedHalfEdges;
    std::unordered_set<Vertex *> verticesToRemove;
    for (const Site &site : mSites)
    {
        HalfEdge *halfEdge = site.cell->outerComponent;
        bool inside = box.contains(halfEdge->origin->point);
        bool outerComponentDirty = !inside;
        HalfEdge *incomingHalfEdge = nullptr; // First half edge coming in the box
        HalfEdge *outgoingHalfEdge = nullptr; // Last half edge going out the box
        Box::Side incomingSide;
        Box::Side outgoingSide;
        do
        {
            std::array<Box::Intersection, 2> intersections;
            int nbIntersections = box.getIntersections(halfEdge->origin->point, halfEdge->destination->point, intersections);
            bool nextInside = box.contains(halfEdge->destination->point);
            HalfEdge *nextHalfEdge = halfEdge->next;
            // The two points are outside the box
            if (!inside && !nextInside)
            {
                // The edge is outside the box
                if (nbIntersections == 0)
                {
                    verticesToRemove.emplace(halfEdge->origin);
                    removeHalfEdge(halfEdge);
                }
                // The edge crosses twice the frontiers of the box
                else if (nbIntersections == 2)
                {
                    verticesToRemove.emplace(halfEdge->origin);
                    if (processedHalfEdges.find(halfEdge->twin) != processedHalfEdges.end())
                    {
                        halfEdge->origin = halfEdge->twin->destination;
                        halfEdge->destination = halfEdge->twin->origin;
                    }
                    else
                    {
                        halfEdge->origin = createVertex(intersections[0].point);
                        halfEdge->destination = createVertex(intersections[1].point);
                    }
                    if (outgoingHalfEdge != nullptr)
                        link(box, outgoingHalfEdge, outgoingSide, halfEdge, intersections[0].side);
                    if (incomingHalfEdge == nullptr)
                    {
                        incomingHalfEdge = halfEdge;
                        incomingSide = intersections[0].side;
                    }
                    outgoingHalfEdge = halfEdge;
                    outgoingSide = intersections[1].side;
                    processedHalfEdges.emplace(halfEdge);
                }
                else
                    success = false;
            }
            // The edge is going outside the box
            else if (inside && !nextInside)
            {
                if (nbIntersections >= 1)
                {
                    if (processedHalfEdges.find(halfEdge->twin) != processedHalfEdges.end())
                        halfEdge->destination = halfEdge->twin->origin;
                    else
                        halfEdge->destination = createVertex(intersections[0].point);
                    outgoingHalfEdge = halfEdge;
                    outgoingSide = intersections[0].side;
                    processedHalfEdges.emplace(halfEdge);
                }
                else
                    success = false;
            }
            // The edge is coming inside the box
            else if (!inside && nextInside)
            {
                if (nbIntersections >= 1)
                {
                    verticesToRemove.emplace(halfEdge->origin);
                    if (processedHalfEdges.find(halfEdge->twin) != processedHalfEdges.end())
                        halfEdge->origin = halfEdge->twin->destination;
                    else
                        halfEdge->origin = createVertex(intersections[0].point);
                    if (outgoingHalfEdge != nullptr)
                        link(box, outgoingHalfEdge, outgoingSide, halfEdge, intersections[0].side);
                    if (incomingHalfEdge == nullptr)
                    {
                        incomingHalfEdge = halfEdge;
                        incomingSide = intersections[0].side;
                    }
                    processedHalfEdges.emplace(halfEdge);
                }
                else
                    success = false;
            }
            halfEdge = nextHalfEdge;
            // Update inside
            inside = nextInside;
        } while (halfEdge != site.cell->outerComponent);
        // Link the last and the first half edges inside the box
        if (outerComponentDirty && incomingHalfEdge != nullptr)
            link(box, outgoingHalfEdge, outgoingSide, incomingHalfEdge, incomingSide);
        // Set outer component
        if (outerComponentDirty)
            site.cell->outerComponent = incomingHalfEdge;
    }
    // Remove vertices
    for (auto &vertex : verticesToRemove)
        removeVertex(vertex);
    // Return the status
    return success;
}

VoronoiDiagram::Vertex *VoronoiDiagram::createVertex(Point point)
{
    mVertices.emplace_back();
    mVertices.back().point = point;
    mVertices.back().it = std::prev(mVertices.end());
    return &mVertices.back();
}

VoronoiDiagram::Vertex *VoronoiDiagram::createCorner(Box box, Box::Side side)
{
    switch (side)
    {
    case Box::Side::LEFT:
        return createVertex(Point(box.left, box.top));
    case Box::Side::BOTTOM:
        return createVertex(Point(box.left, box.bottom));
    case Box::Side::RIGHT:
        return createVertex(Point(box.right, box.bottom));
    case Box::Side::TOP:
        return createVertex(Point(box.right, box.top));
    default:
        return nullptr;
    }
}

VoronoiDiagram::HalfEdge *VoronoiDiagram::createHalfEdge(Cell *cell)
{
    mHalfEdges.emplace_back();
    mHalfEdges.back().incidentFace = cell;
    mHalfEdges.back().it = std::prev(mHalfEdges.end());
    if (cell->outerComponent == nullptr)
        cell->outerComponent = &mHalfEdges.back();
    return &mHalfEdges.back();
}

void VoronoiDiagram::link(Box box, HalfEdge *start, Box::Side startSide, HalfEdge *end, Box::Side endSide)
{
    HalfEdge *halfEdge = start;
    int side = static_cast<int>(startSide);
    while (side != static_cast<int>(endSide))
    {
        side = (side + 1) % 4;
        halfEdge->next = createHalfEdge(start->incidentFace);
        halfEdge->next->prev = halfEdge;
        halfEdge->next->origin = halfEdge->destination;
        halfEdge->next->destination = createCorner(box, static_cast<Box::Side>(side));
        halfEdge = halfEdge->next;
    }
    halfEdge->next = createHalfEdge(start->incidentFace);
    halfEdge->next->prev = halfEdge;
    end->prev = halfEdge->next;
    halfEdge->next->next = end;
    halfEdge->next->origin = halfEdge->destination;
    halfEdge->next->destination = end->origin;
}

void VoronoiDiagram::removeVertex(Vertex *vertex)
{
    mVertices.erase(vertex->it);
}

void VoronoiDiagram::removeHalfEdge(HalfEdge *halfEdge)
{
    mHalfEdges.erase(halfEdge->it);
}


double MAX_BORDER = 1000.0; // handle the far vertices

std::vector<std::vector<Point>> VoronoiDiagram::result()
{
    std::vector<std::vector<Point>> result;
    for (size_t i = 0; i < mSites.size(); i++)
    {
        std::vector<Point> tmp_result;
        const Site *site = &mSites[i];
        Cell *cell = site->cell;
        HalfEdge *halfEdge = cell->outerComponent;
        if (halfEdge == nullptr)
            continue;
        while (halfEdge->prev != nullptr)
        {
            halfEdge = halfEdge->prev;
            if (halfEdge == cell->outerComponent)
                break;
        }
        HalfEdge *start = halfEdge;
        bool far_vertex = false;
        while (halfEdge != nullptr)
        {
            if (halfEdge->origin != nullptr && halfEdge->destination != nullptr)
            {
                if (tmp_result.empty() || (halfEdge->origin->point.getDistance(tmp_result.back()) > 1e-7 && halfEdge->origin->point.getDistance(tmp_result.front()) > 1e-7))
                {
                    tmp_result.push_back(halfEdge->origin->point);
                    if (fabs(halfEdge->origin->point.x) > MAX_BORDER || fabs(halfEdge->origin->point.y) > MAX_BORDER)
                        far_vertex = true;
                }
            }
            halfEdge = halfEdge->next;
            if (halfEdge == start)
                break;
        }

        if (far_vertex)
        {
            std::vector<Point> a = {Point(MAX_BORDER, MAX_BORDER), Point(-MAX_BORDER, MAX_BORDER), Point(-MAX_BORDER, -MAX_BORDER), Point(MAX_BORDER, -MAX_BORDER)};
            tmp_result = convex_clip_polygon(tmp_result, a);
        }

        result.push_back(tmp_result);
    }
    return result;
}

std::vector<std::vector<Point>> VoronoiDiagram::result_with_neighbors()
{
    neighbors.clear();
    std::vector<std::vector<Point>> result;
    for (size_t i = 0; i < mSites.size(); i++)
    {
        std::vector<int> tmp_neighbors;
        std::vector<Point> tmp_result;
        const Site *site = &mSites[i];
        Cell *cell = site->cell;
        HalfEdge *halfEdge = cell->outerComponent;
        if (halfEdge == nullptr)
            continue;
        while (halfEdge->prev != nullptr)
        {
            halfEdge = halfEdge->prev;
            if (halfEdge == cell->outerComponent)
                break;
        }
        HalfEdge *start = halfEdge;
        bool far_vertex = false;
        while (halfEdge != nullptr)
        {
            if (halfEdge->origin != nullptr && halfEdge->destination != nullptr)
            {
                if (tmp_result.empty() || (halfEdge->origin->point.getDistance(tmp_result.back()) > 1e-7 && halfEdge->origin->point.getDistance(tmp_result.front()) > 1e-7))
                {
                    tmp_result.push_back(halfEdge->origin->point);
                    if (fabs(halfEdge->origin->point.x) > MAX_BORDER || fabs(halfEdge->origin->point.y) > MAX_BORDER)
                        far_vertex = true;
                }
            }
            if (halfEdge->twin)
            {
                tmp_neighbors.push_back(halfEdge->twin->incidentFace->site->index);
            }
            else
                tmp_neighbors.push_back(-1);
            halfEdge = halfEdge->next;
            if (halfEdge == start)
                break;
        }
        if (far_vertex)
        {
            std::vector<Point> a = {Point(MAX_BORDER, MAX_BORDER), Point(-MAX_BORDER, MAX_BORDER), Point(-MAX_BORDER, -MAX_BORDER), Point(MAX_BORDER, -MAX_BORDER)};
            tmp_result = convex_clip_polygon(tmp_result, a);
        }
        result.push_back(tmp_result);
        neighbors.push_back(tmp_neighbors);
    }
    return result;
}

std::vector<Point> VoronoiDiagram::result(size_t i)
{
    std::vector<Point> result;
    const Site *site = &mSites[i];
    Cell *cell = site->cell;
    HalfEdge *halfEdge = cell->outerComponent;
    HalfEdge *start = halfEdge;
    bool far_vertex = false;
    while (halfEdge != nullptr)
    {
        if (halfEdge->origin != nullptr && halfEdge->destination != nullptr)
        {
            result.push_back(halfEdge->origin->point);
            if (fabs(halfEdge->origin->point.x) > MAX_BORDER || fabs(halfEdge->origin->point.y) > MAX_BORDER)
                far_vertex = true;
        }
        halfEdge = halfEdge->next;
        if (halfEdge == start)
            break;
    }
    if (far_vertex)
    {
        std::vector<Point> a = {Point(MAX_BORDER, MAX_BORDER), Point(-MAX_BORDER, MAX_BORDER), Point(-MAX_BORDER, -MAX_BORDER), Point(MAX_BORDER, -MAX_BORDER)};
        result = convex_clip_polygon(result, a);
    }
    return result;
}

// Event
Event::Event(VoronoiDiagram::Site *site) : type(Type::SITE), x(site->point.x), y(site->point.y), index(-1), site(site)
{
}

Event::Event(double y, Point point, Arc *arc) : type(Type::CIRCLE), x(point.x), y(y), index(-1), point(point), arc(arc)
{
}

bool operator<(const Event &lhs, const Event &rhs)
{
    return lhs.y < rhs.y || (fabs(lhs.y - rhs.y) <= 1e-10 && lhs.x < rhs.x);
}

// Beachline
Beachline::Beachline() : mNil(new Arc), mRoot(mNil)
{
    mNil->color = Arc::Color::BLACK;
}

Beachline::~Beachline()
{
    free(mRoot);
    delete mNil;
}

Beachline::Beachline(Beachline &&other)
{
    mNil = other.mNil;
    mRoot = other.mRoot;
    other.mNil = nullptr;
    other.mRoot = nullptr;
}

Beachline &Beachline::operator=(Beachline &&other)
{
    free(mRoot);
    delete mNil;
    mNil = other.mNil;
    mRoot = other.mRoot;
    other.mNil = nullptr;
    other.mRoot = nullptr;
    return *this;
}

Arc *Beachline::createArc(VoronoiDiagram::Site *site, Arc::Side side)
{
    return new Arc{mNil, mNil, mNil, site, nullptr, nullptr, nullptr, mNil, mNil, Arc::Color::RED, side};
}

bool Beachline::isEmpty() const
{
    return isNil(mRoot);
}

bool Beachline::isNil(const Arc *x) const
{
    return x == mNil;
}

void Beachline::setRoot(Arc *x)
{
    mRoot = x;
    mRoot->color = Arc::Color::BLACK;
}

Arc *Beachline::getLeftmostArc() const
{
    Arc *x = mRoot;
    while (!isNil(x->prev))
        x = x->prev;
    return x;
}

Arc *Beachline::locateArcAbove(const Point &point, double l) const
{
    Arc *node = mRoot;
    bool found = false;
    while (!found)
    {
        double breakpointLeft = -std::numeric_limits<double>::infinity();
        double breakpointRight = std::numeric_limits<double>::infinity();
        if (!isNil(node->prev))
            breakpointLeft = computeBreakpoint(node->prev->site->point, node->site->point, l, node->prev->side);
        if (!isNil(node->next))
            breakpointRight = computeBreakpoint(node->site->point, node->next->site->point, l, node->next->side);
        if (point.x < breakpointLeft)
            node = node->left;
        else if (point.x > breakpointRight)
            node = node->right;
        else
            found = true;
    }
    return node;
}

void Beachline::insertBefore(Arc *x, Arc *y)
{
    // Find the right place
    if (isNil(x->left))
    {
        x->left = y;
        y->parent = x;
    }
    else
    {
        x->prev->right = y;
        y->parent = x->prev;
    }
    // Set the pointers
    y->prev = x->prev;
    if (!isNil(y->prev))
        y->prev->next = y;
    y->next = x;
    x->prev = y;
    // Balance the tree
    insertFixup(y);
}

void Beachline::insertAfter(Arc *x, Arc *y)
{
    // Find the right place
    if (isNil(x->right))
    {
        x->right = y;
        y->parent = x;
    }
    else
    {
        x->next->left = y;
        y->parent = x->next;
    }
    // Set the pointers
    y->next = x->next;
    if (!isNil(y->next))
        y->next->prev = y;
    y->prev = x;
    x->next = y;
    // Balance the tree
    insertFixup(y);
}

void Beachline::replace(Arc *x, Arc *y)
{
    transplant(x, y);
    y->left = x->left;
    y->right = x->right;
    if (!isNil(y->left))
        y->left->parent = y;
    if (!isNil(y->right))
        y->right->parent = y;
    y->prev = x->prev;
    y->next = x->next;
    if (!isNil(y->prev))
        y->prev->next = y;
    if (!isNil(y->next))
        y->next->prev = y;
    y->color = x->color;
}

void Beachline::remove(Arc *z)
{
    Arc *y = z;
    Arc::Color yOriginalColor = y->color;
    Arc *x = mNil;
    if (isNil(z->left))
    {
        x = z->right;
        transplant(z, z->right);
    }
    else if (isNil(z->right))
    {
        x = z->left;
        transplant(z, z->left);
    }
    else
    {
        y = minimum(z->right);
        yOriginalColor = y->color;
        x = y->right;
        if (y->parent == z)
            x->parent = y; // Because x could be Nil
        else
        {
            transplant(y, y->right);
            y->right = z->right;
            y->right->parent = y;
        }
        transplant(z, y);
        y->left = z->left;
        y->left->parent = y;
        y->color = z->color;
    }
    if (yOriginalColor == Arc::Color::BLACK)
        removeFixup(x);
    // Update next and prev
    if (!isNil(z->prev))
        z->prev->next = z->next;
    if (!isNil(z->next))
        z->next->prev = z->prev;
}

Arc *Beachline::minimum(Arc *x) const
{
    while (!isNil(x->left))
        x = x->left;
    return x;
}

void Beachline::transplant(Arc *u, Arc *v)
{
    if (isNil(u->parent))
        mRoot = v;
    else if (u == u->parent->left)
        u->parent->left = v;
    else
        u->parent->right = v;
    v->parent = u->parent;
}

void Beachline::insertFixup(Arc *z)
{
    while (z->parent->color == Arc::Color::RED)
    {
        if (z->parent == z->parent->parent->left)
        {
            Arc *y = z->parent->parent->right;
            // Case 1
            if (y->color == Arc::Color::RED)
            {
                z->parent->color = Arc::Color::BLACK;
                y->color = Arc::Color::BLACK;
                z->parent->parent->color = Arc::Color::RED;
                z = z->parent->parent;
            }
            else
            {
                // Case 2
                if (z == z->parent->right)
                {
                    z = z->parent;
                    leftRotate(z);
                }
                // Case 3
                z->parent->color = Arc::Color::BLACK;
                z->parent->parent->color = Arc::Color::RED;
                rightRotate(z->parent->parent);
            }
        }
        else
        {
            Arc *y = z->parent->parent->left;
            // Case 1
            if (y->color == Arc::Color::RED)
            {
                z->parent->color = Arc::Color::BLACK;
                y->color = Arc::Color::BLACK;
                z->parent->parent->color = Arc::Color::RED;
                z = z->parent->parent;
            }
            else
            {
                // Case 2
                if (z == z->parent->left)
                {
                    z = z->parent;
                    rightRotate(z);
                }
                // Case 3
                z->parent->color = Arc::Color::BLACK;
                z->parent->parent->color = Arc::Color::RED;
                leftRotate(z->parent->parent);
            }
        }
    }
    mRoot->color = Arc::Color::BLACK;
}

void Beachline::removeFixup(Arc *x)
{
    while (x != mRoot && x->color == Arc::Color::BLACK)
    {
        if (x == x->parent->left)
        {
            Arc *w = x->parent->right;
            // Case 1
            if (w->color == Arc::Color::RED)
            {
                w->color = Arc::Color::BLACK;
                x->parent->color = Arc::Color::RED;
                leftRotate(x->parent);
                w = x->parent->right;
            }
            // Case 2
            if (w->left->color == Arc::Color::BLACK && w->right->color == Arc::Color::BLACK)
            {
                w->color = Arc::Color::RED;
                x = x->parent;
            }
            else
            {
                // Case 3
                if (w->right->color == Arc::Color::BLACK)
                {
                    w->left->color = Arc::Color::BLACK;
                    w->color = Arc::Color::RED;
                    rightRotate(w);
                    w = x->parent->right;
                }
                // Case 4
                w->color = x->parent->color;
                x->parent->color = Arc::Color::BLACK;
                w->right->color = Arc::Color::BLACK;
                leftRotate(x->parent);
                x = mRoot;
            }
        }
        else
        {
            Arc *w = x->parent->left;
            // Case 1
            if (w->color == Arc::Color::RED)
            {
                w->color = Arc::Color::BLACK;
                x->parent->color = Arc::Color::RED;
                rightRotate(x->parent);
                w = x->parent->left;
            }
            // Case 2
            if (w->left->color == Arc::Color::BLACK && w->right->color == Arc::Color::BLACK)
            {
                w->color = Arc::Color::RED;
                x = x->parent;
            }
            else
            {
                // Case 3
                if (w->left->color == Arc::Color::BLACK)
                {
                    w->right->color = Arc::Color::BLACK;
                    w->color = Arc::Color::RED;
                    leftRotate(w);
                    w = x->parent->left;
                }
                // Case 4
                w->color = x->parent->color;
                x->parent->color = Arc::Color::BLACK;
                w->left->color = Arc::Color::BLACK;
                rightRotate(x->parent);
                x = mRoot;
            }
        }
    }
    x->color = Arc::Color::BLACK;
}

void Beachline::leftRotate(Arc *x)
{
    Arc *y = x->right;
    x->right = y->left;
    if (!isNil(y->left))
        y->left->parent = x;
    y->parent = x->parent;
    if (isNil(x->parent))
        mRoot = y;
    else if (x->parent->left == x)
        x->parent->left = y;
    else
        x->parent->right = y;
    y->left = x;
    x->parent = y;
}

void Beachline::rightRotate(Arc *y)
{
    Arc *x = y->left;
    y->left = x->right;
    if (!isNil(x->right))
        x->right->parent = y;
    x->parent = y->parent;
    if (isNil(y->parent))
        mRoot = x;
    else if (y->parent->left == y)
        y->parent->left = x;
    else
        y->parent->right = x;
    x->right = y;
    y->parent = x;
}

double Beachline::computeBreakpoint(const Point &point1, const Point &point2, double l, Arc::Side side) const
{
    double x1 = point1.x;
    double y1 = point1.y;
    double x2 = point2.x;
    double y2 = point2.y;
    if (fabs(y1 - y2) <= 1e-10)
    {
        // The breakpoint is between the two points
        if (x1 < x2)
            return (x1 + x2) / 2;
        // The breakpoint is at infinity
        else
            return side == Arc::Side::Left ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
    }
    if (fabs(y1 - l) <= 1e-10)
        return x1;
    if (fabs(y2 - l) <= 1e-10)
        return x2;
    double d1 = 1.0 / (2.0 * (y1 - l));
    double d2 = 1.0 / (2.0 * (y2 - l));
    double a = d1 - d2;
    double b = 2.0 * (x2 * d2 - x1 * d1);
    double c = (y1 * y1 + x1 * x1 - l * l) * d1 - (y2 * y2 + x2 * x2 - l * l) * d2;
    double delta = b * b - 4.0 * a * c;
    return (-b + std::sqrt(delta)) / (2.0 * a);
}

void Beachline::free(Arc *x)
{
    if (isNil(x))
        return;
    else
    {
        free(x->left);
        free(x->right);
        delete x;
    }
}

//  FortuneAlgorithm
FortuneAlgorithm::FortuneAlgorithm(std::vector<Point> points) : mDiagram(std::move(points))
{
}

FortuneAlgorithm::~FortuneAlgorithm() = default;

void FortuneAlgorithm::construct()
{
    // Initialize event queue
    for (std::size_t i = 0; i < mDiagram.getNbSites(); i++)
        mEvents.push(std::make_unique<Event>(mDiagram.getSite(i)));

    // Process events
    while (!mEvents.isEmpty())
    {
        std::unique_ptr<Event> event = mEvents.pop();
        mBeachlineY = event->y;
        if (event->type == Event::Type::SITE)
            handleSiteEvent(event.get());
        else
            handleCircleEvent(event.get());
    }
}

VoronoiDiagram FortuneAlgorithm::getDiagram()
{
    return std::move(mDiagram);
}

void FortuneAlgorithm::handleSiteEvent(Event *event)
{
    VoronoiDiagram::Site *site = event->site;
    // 1. Check if the bachline is empty
    if (mBeachline.isEmpty())
    {
        mBeachline.setRoot(mBeachline.createArc(site, Arc::Side::Left));
        return;
    }
    // 2. Look for the arc above the site
    Arc *arcToBreak = mBeachline.locateArcAbove(site->point, mBeachlineY);
    deleteEvent(arcToBreak);
    // 3. Replace this arc by the new arcs
    Arc *middleArc = breakArc(arcToBreak, site);

    Arc *leftArc = middleArc->prev;
    Arc *rightArc = middleArc->next;

    // 4. Add an edge in the diagram
    addEdge(leftArc, middleArc);
    middleArc->rightHalfEdge = middleArc->leftHalfEdge;
    rightArc->leftHalfEdge = leftArc->rightHalfEdge;
    // 5. Check circle events
    // Left triplet
    if (!mBeachline.isNil(leftArc->prev))
        addEvent(leftArc->prev, leftArc, middleArc);
    // Right triplet
    if (!mBeachline.isNil(rightArc->next))
        addEvent(middleArc, rightArc, rightArc->next);
}

void FortuneAlgorithm::handleCircleEvent(Event *event)
{
    Point point = event->point;
    Arc *arc = event->arc;
    // 1. Add vertex
    VoronoiDiagram::Vertex *vertex = mDiagram.createVertex(point);
    // 2. Delete all the events with this arc
    Arc *leftArc = arc->prev;
    Arc *rightArc = arc->next;
    deleteEvent(leftArc);
    deleteEvent(rightArc);
    // 3. Update the beachline and the diagram
    removeArc(arc, vertex);
    // 4. Add new circle events
    // Left triplet
    if (!mBeachline.isNil(leftArc->prev))
        addEvent(leftArc->prev, leftArc, rightArc);
    // Right triplet
    if (!mBeachline.isNil(rightArc->next))
        addEvent(leftArc, rightArc, rightArc->next);
}

Arc *FortuneAlgorithm::breakArc(Arc *arc, VoronoiDiagram::Site *site)
{
    // Create the new subtree
    Arc *middleArc = mBeachline.createArc(site, Arc::Side::Left);
    Arc *leftArc = mBeachline.createArc(arc->site, Arc::Side::Left);
    leftArc->leftHalfEdge = arc->leftHalfEdge;
    Arc *rightArc = mBeachline.createArc(arc->site, Arc::Side::Right);
    rightArc->rightHalfEdge = arc->rightHalfEdge;
    // Insert the subtree in the beachline
    mBeachline.replace(arc, middleArc);
    mBeachline.insertBefore(middleArc, leftArc);
    mBeachline.insertAfter(middleArc, rightArc);
    // Delete old arc
    delete arc;
    // Return the middle arc
    return middleArc;
}

void FortuneAlgorithm::removeArc(Arc *arc, VoronoiDiagram::Vertex *vertex)
{
    // End edges
    setDestination(arc->prev, arc, vertex);
    setDestination(arc, arc->next, vertex);
    // Join the edges of the middle arc
    arc->leftHalfEdge->next = arc->rightHalfEdge;
    arc->rightHalfEdge->prev = arc->leftHalfEdge;
    // Update beachline
    mBeachline.remove(arc);
    // Create a new edge
    VoronoiDiagram::HalfEdge *prevHalfEdge = arc->prev->rightHalfEdge;
    VoronoiDiagram::HalfEdge *nextHalfEdge = arc->next->leftHalfEdge;
    addEdge(arc->prev, arc->next);
    setOrigin(arc->prev, arc->next, vertex);
    setPrevHalfEdge(arc->prev->rightHalfEdge, prevHalfEdge);
    setPrevHalfEdge(nextHalfEdge, arc->next->leftHalfEdge);
    // Delete node
    delete arc;
}

bool FortuneAlgorithm::isMovingRight(const Arc *left, const Arc *right) const
{
    return left->site->point.y < right->site->point.y - 1e-10;
}

double FortuneAlgorithm::getInitialX(const Arc *left, const Arc *right, bool movingRight) const
{
    return movingRight ? left->site->point.x : right->site->point.x;
}

void FortuneAlgorithm::addEdge(Arc *left, Arc *right)
{
    // Create two new half edges
    left->rightHalfEdge = mDiagram.createHalfEdge(left->site->cell);
    right->leftHalfEdge = mDiagram.createHalfEdge(right->site->cell);
    // Set the two half edges twins
    left->rightHalfEdge->twin = right->leftHalfEdge;
    right->leftHalfEdge->twin = left->rightHalfEdge;
}

void FortuneAlgorithm::setOrigin(Arc *left, Arc *right, VoronoiDiagram::Vertex *vertex)
{
    left->rightHalfEdge->destination = vertex;
    right->leftHalfEdge->origin = vertex;
}

void FortuneAlgorithm::setDestination(Arc *left, Arc *right, VoronoiDiagram::Vertex *vertex)
{
    left->rightHalfEdge->origin = vertex;
    right->leftHalfEdge->destination = vertex;
}

void FortuneAlgorithm::setPrevHalfEdge(VoronoiDiagram::HalfEdge *prev, VoronoiDiagram::HalfEdge *next)
{
    prev->next = next;
    next->prev = prev;
}

void FortuneAlgorithm::addEvent(Arc *left, Arc *middle, Arc *right)
{
    double y;
    Point convergencePoint = computeConvergencePoint(left->site->point, middle->site->point, right->site->point, y);

    // Check that the bottom of the circle is below the beachline
    if (!(y <= mBeachlineY + 1e-10))
        return;
    // Check that the event is valid by checking that the middle arc will indeed shrink
    // To do that we check in which direction the breakpoint are moving and that the
    // convergence point is somewhere the breakpoints are moving to
    bool leftBreakpointMovingRight = isMovingRight(left, middle);
    bool rightBreakpointMovingRight = isMovingRight(middle, right);
    if (fabs(left->site->point.y - middle->site->point.y) <= 1e-10)
    {
        leftBreakpointMovingRight = true;
    }

    if (fabs(middle->site->point.y - right->site->point.y) <= 1e-10)
    {
        rightBreakpointMovingRight = false;
    }

    double leftInitialX = getInitialX(left, middle, leftBreakpointMovingRight);
    double rightInitialX = getInitialX(middle, right, rightBreakpointMovingRight);

    bool isValid =
        !(!leftBreakpointMovingRight && rightBreakpointMovingRight) &&
        ((leftBreakpointMovingRight && (leftInitialX <= convergencePoint.x + 1e-6)) ||
         (!leftBreakpointMovingRight && (leftInitialX >= convergencePoint.x - 1e-6))) &&
        ((rightBreakpointMovingRight && (rightInitialX <= convergencePoint.x + 1e-6)) ||
         (!rightBreakpointMovingRight && (rightInitialX >= convergencePoint.x - 1e-6)));
    if (isValid)
    {
        std::unique_ptr<Event> event = std::make_unique<Event>(y, convergencePoint, middle);
        middle->event = event.get();
        mEvents.push(std::move(event));
    }
}

void FortuneAlgorithm::deleteEvent(Arc *arc)
{
    if (arc->event != nullptr)
    {
        mEvents.remove(arc->event->index);
        arc->event = nullptr;
    }
}

Point FortuneAlgorithm::computeConvergencePoint(const Point &point1, const Point &point2, const Point &point3, double &y) const
{
    Point v1 = (point1 - point2).getOrthogonal();
    Point v2 = (point2 - point3).getOrthogonal();
    Point delta = 0.5 * (point3 - point1);
    if (fabs(v1.getDet(v2)) < 1e-8)
    {
        // Infinity means that the event will never be added to the event queue
        y = std::numeric_limits<double>::infinity();
        return Point(0.0, 0.0);
    }
    double t = delta.getDet(v2) / v1.getDet(v2);
    Point center = 0.5 * (point1 + point2) + t * v1;
    double r = center.getDistance(point1);
    y = center.y - r;
    return center;
}

bool FortuneAlgorithm::bound(Box box)
{
    bool success = true;
    // Make sure the bounding box contains all the vertices
    for (const auto &vertex : mDiagram.getVertices()) // Much faster when using vector<unique_ptr<Vertex*>, maybe we can test vertices in border cells to speed up
    {
        box.left = std::min(vertex.point.x, box.left);
        box.bottom = std::min(vertex.point.y, box.bottom);
        box.right = std::max(vertex.point.x, box.right);
        box.top = std::max(vertex.point.y, box.top);
    }
    // Retrieve all non bounded half edges from the beachline
    std::list<LinkedVertex> linkedVertices;
    std::unordered_map<std::size_t, std::array<LinkedVertex *, 8>> vertices(mDiagram.getNbSites());
    if (!mBeachline.isEmpty())
    {
        Arc *leftArc = mBeachline.getLeftmostArc();
        Arc *rightArc = leftArc->next;
        while (!mBeachline.isNil(rightArc))
        {
            // Bound the edge
            Point direction = (leftArc->site->point - rightArc->site->point).getOrthogonal();
            Point origin = (leftArc->site->point + rightArc->site->point) * 0.5;
            // Line-box intersection
            Box::Intersection intersection = box.getFirstIntersection(origin, direction);
            // Create a new vertex and ends the half edges
            VoronoiDiagram::Vertex *vertex = mDiagram.createVertex(intersection.point);
            setDestination(leftArc, rightArc, vertex);
            // Initialize pointers
            if (vertices.find(leftArc->site->index) == vertices.end())
                vertices[leftArc->site->index].fill(nullptr);
            if (vertices.find(rightArc->site->index) == vertices.end())
                vertices[rightArc->site->index].fill(nullptr);

            // Check that the vertices are not already assigned
            success = vertices[leftArc->site->index][2 * static_cast<int>(intersection.side) + 1] == nullptr && success;
            success = vertices[rightArc->site->index][2 * static_cast<int>(intersection.side)] == nullptr && success;
            // Store the vertex on the boundaries
            linkedVertices.emplace_back(LinkedVertex{nullptr, vertex, leftArc->rightHalfEdge});
            vertices[leftArc->site->index][2 * static_cast<int>(intersection.side) + 1] = &linkedVertices.back();
            linkedVertices.emplace_back(LinkedVertex{rightArc->leftHalfEdge, vertex, nullptr});
            vertices[rightArc->site->index][2 * static_cast<int>(intersection.side)] = &linkedVertices.back();
            // Next edge
            leftArc = rightArc;
            rightArc = rightArc->next;
        }
    }
    // Add corners
    for (auto &kv : vertices)
    {
        auto &cellVertices = kv.second;
        // We check twice the first side to be sure that all necessary corners are added
        for (std::size_t i = 0; i < 5; i++)
        {
            std::size_t side = i % 4;
            std::size_t nextSide = (side + 1) % 4;
            // Add first corner
            if (cellVertices[2 * side] == nullptr && cellVertices[2 * side + 1] != nullptr)
            {
                std::size_t prevSide = (side + 3) % 4;
                VoronoiDiagram::Vertex *corner = mDiagram.createCorner(box, static_cast<Box::Side>(side));
                linkedVertices.emplace_back(LinkedVertex{nullptr, corner, nullptr});
                // Check that we are not erasing an already assigned vertex
                success = cellVertices[2 * prevSide + 1] == nullptr && success;
                // Store the vertex on the boundary
                cellVertices[2 * prevSide + 1] = &linkedVertices.back();
                cellVertices[2 * side] = &linkedVertices.back();
            }
            // Add second corner
            else if (cellVertices[2 * side] != nullptr && cellVertices[2 * side + 1] == nullptr)
            {
                VoronoiDiagram::Vertex *corner = mDiagram.createCorner(box, static_cast<Box::Side>(nextSide));
                linkedVertices.emplace_back(LinkedVertex{nullptr, corner, nullptr});
                // Check that we are not erasing an already assigned vertex
                success = cellVertices[2 * nextSide] == nullptr && success;
                // Store the vertex on the boundary
                cellVertices[2 * side + 1] = &linkedVertices.back();
                cellVertices[2 * nextSide] = &linkedVertices.back();
            }
        }
    }
    // Join the half edges
    for (auto &kv : vertices)
    {
        std::size_t i = kv.first;
        auto &cellVertices = kv.second;
        for (std::size_t side = 0; side < 4; side++)
        {
            if (cellVertices[2 * side] != nullptr)
            {
                // Link vertices
                VoronoiDiagram::HalfEdge *halfEdge = mDiagram.createHalfEdge(mDiagram.getCell(i));
                halfEdge->origin = cellVertices[2 * side]->vertex;
                halfEdge->destination = cellVertices[2 * side + 1]->vertex;
                cellVertices[2 * side]->nextHalfEdge = halfEdge;
                halfEdge->prev = cellVertices[2 * side]->prevHalfEdge;
                if (cellVertices[2 * side]->prevHalfEdge != nullptr)
                    cellVertices[2 * side]->prevHalfEdge->next = halfEdge;
                cellVertices[2 * side + 1]->prevHalfEdge = halfEdge;
                halfEdge->next = cellVertices[2 * side + 1]->nextHalfEdge;
                if (cellVertices[2 * side + 1]->nextHalfEdge != nullptr)
                    cellVertices[2 * side + 1]->nextHalfEdge->prev = halfEdge;
            }
        }
    }
    return success; // TO DO: detect errors
}
