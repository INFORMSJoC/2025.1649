#pragma once

#include <cstdlib>
#include <stdint.h> //#541
#include <iostream>
#include <queue>
#include <vector>
#include <functional>
#include <numeric>
#include <memory>

#include "clipper.core.h"

namespace Clipper2Lib
{

	struct Scanline;
	struct IntersectNode;
	struct Active;
	struct Vertex;
	struct LocalMinima;
	struct OutRec;
	struct HorzSegment;

	// Note: all clipping operations except for Difference are commutative.
	enum class ClipType
	{
		None,
		Intersection,
		Union,
		Difference,
		Xor
	};

	enum class PathType
	{
		Subject,
		Clip
	};

	enum class JoinWith
	{
		None,
		Left,
		Right
	};

	enum class VertexFlags : uint32_t
	{
		None = 0,
		OpenStart = 1,
		OpenEnd = 2,
		LocalMax = 4,
		LocalMin = 8
	};

	constexpr enum VertexFlags operator&(enum VertexFlags a, enum VertexFlags b)
	{
		return (enum VertexFlags)(uint32_t(a) & uint32_t(b));
	}

	constexpr enum VertexFlags operator|(enum VertexFlags a, enum VertexFlags b)
	{
		return (enum VertexFlags)(uint32_t(a) | uint32_t(b));
	}

	struct Vertex
	{
		Point pt;
		Vertex *next = nullptr;
		Vertex *prev = nullptr;
		VertexFlags flags = VertexFlags::None;
	};

	struct OutPt
	{
		Point pt;
		OutPt *next = nullptr;
		OutPt *prev = nullptr;
		OutRec *outrec;
		HorzSegment *horz = nullptr;

		OutPt(const Point &pt_, OutRec *outrec_) : pt(pt_), outrec(outrec_)
		{
			next = this;
			prev = this;
		}
	};

	class PolyPaths;
	class PolyPath;
	using PolyTree = PolyPath;
	struct OutRec;
	typedef std::vector<OutRec *> OutRecList;

	// OutRec: contains a path in the clipping solution. Edges in the AEL will
	// have OutRec pointers assigned when they form part of the clipping solution.
	struct OutRec
	{
		size_t idx = 0;
		OutRec *owner = nullptr;
		Active *front_edge = nullptr;
		Active *back_edge = nullptr;
		OutPt *pts = nullptr;
		PolyPaths *polypath = nullptr;
		OutRecList *splits = nullptr;
		OutRec *recursive_split = nullptr;
		Rect bounds = {};
		Path path;
		bool is_open = false;

		~OutRec()
		{
			if (splits)
				delete splits;
			// nb: don't delete the split pointers
			// as these are owned by ClipperBase's outrec_list_
		};
	};

	///////////////////////////////////////////////////////////////////
	// Important: UP and DOWN here are premised on Y-axis positive down
	// displays, which is the orientation used in Clipper's development.
	///////////////////////////////////////////////////////////////////

	struct Active
	{
		Point bot;
		Point top;
		int64_t curr_x = 0; // current (updated at every new scanline)
		double dx = 0.0;
		int wind_dx = 1; // 1 or -1 depending on winding direction
		int wind_cnt = 0;
		int wind_cnt2 = 0; // winding count of the opposite polytype
		OutRec *outrec = nullptr;
		// AEL: 'active edge list' (Vatti's AET - active edge table)
		//      a linked list of all edges (from left to right) that are present
		//      (or 'active') within the current scanbeam (a horizontal 'beam' that
		//      sweeps from bottom to top over the paths in the clipping operation).
		Active *prev_in_ael = nullptr;
		Active *next_in_ael = nullptr;
		// SEL: 'sorted edge list' (Vatti's ST - sorted table)
		//      linked list used when sorting edges into their new positions at the
		//      top of scanbeams, but also (re)used to process horizontals.
		Active *prev_in_sel = nullptr;
		Active *next_in_sel = nullptr;
		Active *jump = nullptr;
		Vertex *vertex_top = nullptr;
		LocalMinima *local_min = nullptr; // the bottom of an edge 'bound' (also Vatti)
		bool is_left_bound = false;
		JoinWith join_with = JoinWith::None;
	};

	struct LocalMinima
	{
		Vertex *vertex;
		PathType polytype;
		bool is_open;
		LocalMinima(Vertex *v, PathType pt, bool open) : vertex(v), polytype(pt), is_open(open) {}
	};

	struct IntersectNode
	{
		Point pt;
		Active *edge1;
		Active *edge2;
		IntersectNode() : pt(Point(0, 0)), edge1(NULL), edge2(NULL) {}
		IntersectNode(Active *e1, Active *e2, Point &pt_) : pt(pt_), edge1(e1), edge2(e2) {}
	};

	struct HorzSegment
	{
		OutPt *left_op;
		OutPt *right_op = nullptr;
		bool left_to_right = true;
		HorzSegment() : left_op(nullptr) {}
		explicit HorzSegment(OutPt *op) : left_op(op) {}
	};

	struct HorzJoin
	{
		OutPt *op1 = nullptr;
		OutPt *op2 = nullptr;
		HorzJoin(){};
		explicit HorzJoin(OutPt *ltr, OutPt *rtl) : op1(ltr), op2(rtl) {}
	};

	typedef std::vector<HorzSegment> HorzSegmentList;
	typedef std::unique_ptr<LocalMinima> LocalMinima_ptr;
	typedef std::vector<LocalMinima_ptr> LocalMinimaList;
	typedef std::vector<IntersectNode> IntersectNodeList;

	// ReuseableDataContainer ------------------------------------------------

	class ReuseableDataContainer
	{
	private:
		friend class ClipperBase;
		LocalMinimaList minima_list_;
		std::vector<Vertex *> vertex_lists_;
		void AddLocMin(Vertex &vert, PathType polytype, bool is_open);

	public:
		virtual ~ReuseableDataContainer();
		void Clear();
		void AddPaths(const Paths &paths, PathType polytype, bool is_open);
	};

	// ClipperBase -------------------------------------------------------------

	class ClipperBase
	{
	private:
		ClipType cliptype_ = ClipType::None;
		FillRule fillrule_ = FillRule::EvenOdd;
		FillRule fillpos = FillRule::Positive;
		int64_t bot_y_ = 0;
		bool minima_list_sorted_ = false;
		bool using_polytree_ = false;
		Active *actives_ = nullptr;
		Active *sel_ = nullptr;
		LocalMinimaList minima_list_; // pointers in case of memory reallocs
		LocalMinimaList::iterator current_locmin_iter_;
		std::vector<Vertex *> vertex_lists_;
		std::priority_queue<int64_t> scanline_list_;
		IntersectNodeList intersect_nodes_;
		HorzSegmentList horz_seg_list_;
		std::vector<HorzJoin> horz_join_list_;
		void Reset();
		inline void InsertScanline(int64_t y);
		inline bool PopScanline(int64_t &y);
		inline bool PopLocalMinima(int64_t y, LocalMinima *&local_minima);
		void DisposeAllOutRecs();
		void DisposeVerticesAndLocalMinima();
		void DeleteEdges(Active *&e);
		inline void AddLocMin(Vertex &vert, PathType polytype, bool is_open);
		bool IsContributingClosed(const Active &e) const;
		inline bool IsContributingOpen(const Active &e) const;
		void SetWindCountForClosedPathEdge(Active &edge);
		void SetWindCountForOpenPathEdge(Active &e);
		void InsertLocalMinimaIntoAEL(int64_t bot_y);
		void InsertLeftEdge(Active &e);
		inline void PushHorz(Active &e);
		inline bool PopHorz(Active *&e);
		inline OutPt *StartOpenPath(Active &e, const Point &pt);
		inline void UpdateEdgeIntoAEL(Active *e);
		OutPt *IntersectEdges(Active &e1, Active &e2, const Point &pt);
		inline void DeleteFromAEL(Active &e);
		inline void AdjustCurrXAndCopyToSEL(const int64_t top_y);
		void DoIntersections(const int64_t top_y);
		void AddNewIntersectNode(Active &e1, Active &e2, const int64_t top_y);
		bool BuildIntersectList(const int64_t top_y);
		void ProcessIntersectList();
		void SwapPositionsInAEL(Active &edge1, Active &edge2);
		OutRec *NewOutRec();
		OutPt *AddOutPt(const Active &e, const Point &pt);
		OutPt *AddLocalMinPoly(Active &e1, Active &e2, const Point &pt, bool is_new = false);
		OutPt *AddLocalMaxPoly(Active &e1, Active &e2, const Point &pt);
		void DoHorizontal(Active &horz);
		bool ResetHorzDirection(const Active &horz, const Vertex *max_vertex, int64_t &horz_left, int64_t &horz_right);
		void DoTopOfScanbeam(const int64_t top_y);
		Active *DoMaxima(Active &e);
		void JoinOutrecPaths(Active &e1, Active &e2);
		void FixSelfIntersects(OutRec *outrec);
		void DoSplitOp(OutRec *outRec, OutPt *splitOp);

		inline void AddTrialHorzJoin(OutPt *op);
		void ConvertHorzSegsToJoins();
		void ProcessHorzJoins();

		void Split(Active &e, const Point &pt);
		inline void CheckJoinLeft(Active &e, const Point &pt, bool check_curr_x = false);
		inline void CheckJoinRight(Active &e, const Point &pt, bool check_curr_x = false);

	protected:
		int error_code_ = 0;
		bool has_open_paths_ = false;
		bool succeeded_ = true;
		OutRecList outrec_list_; // pointers in case list memory reallocated
		bool ExecuteInternal(ClipType ct, FillRule ft, bool use_polytrees);
		void CleanCollinear(OutRec *outrec);
		bool CheckBounds(OutRec *outrec);
		bool CheckSplitOwner(OutRec *outrec, OutRecList *splits);
		void RecursiveCheckOwners(OutRec *outrec, PolyPaths *polypath);
		void CleanUp(); // unlike Clear, CleanUp preserves added paths
		void AddPath(const Path &path, PathType polytype, bool is_open);
		void AddPaths(const Paths &paths, PathType polytype, bool is_open);

	public:
		virtual ~ClipperBase();
		int ErrorCode() { return error_code_; };
		bool PreserveCollinear = true;
		bool ReverseSolution = false;
		void Clear();
		void AddReuseableData(const ReuseableDataContainer &reuseable_data);
	};

	// PolyPath / PolyTree --------------------------------------------------------

	// PolyTree: is intended as a READ-ONLY data structure for CLOSED paths returned
	// by clipping operations. While this structure is more complex than the
	// alternative Paths structure, it does preserve path 'ownership' - ie those
	// paths that contain (or own) other paths. This will be useful to some users.

	class PolyPaths
	{
	protected:
		PolyPaths *parent_;

	public:
		PolyPaths(PolyPaths *parent = nullptr) : parent_(parent) {}
		virtual ~PolyPaths(){};
		// https://en.cppreference.com/w/cpp/language/rule_of_three
		PolyPaths(const PolyPaths &) = delete;
		PolyPaths &operator=(const PolyPaths &) = delete;

		unsigned Level() const
		{
			unsigned result = 0;
			const PolyPaths *p = parent_;
			while (p)
			{
				++result;
				p = p->parent_;
			}
			return result;
		}

		virtual PolyPaths *AddChild(const Path &path) = 0;

		virtual void Clear() = 0;
		virtual size_t Count() const { return 0; }

		const PolyPaths *Parent() const { return parent_; }

		bool IsHole() const
		{
			unsigned lvl = Level();
			// Even levels except level 0
			return lvl && !(lvl & 1);
		}
	};

	typedef typename std::vector<std::unique_ptr<PolyPath>> PolyPathDList;

	class PolyPath : public PolyPaths
	{
	private:
		PolyPathDList childs_;
		double scale_;
		Path polygon_;

	public:
		explicit PolyPath(PolyPath *parent = nullptr) : PolyPaths(parent)
		{
			scale_ = parent ? parent->scale_ : 1.0;
		}

		~PolyPath()
		{
			childs_.resize(0);
		}

		const PolyPath *operator[](size_t index) const
		{
			return childs_[index].get();
		}

		const PolyPath *Child(size_t index) const
		{
			return childs_[index].get();
		}

		PolyPathDList::const_iterator begin() const { return childs_.cbegin(); }
		PolyPathDList::const_iterator end() const { return childs_.cend(); }

		void SetScale(double value) { scale_ = value; }
		double Scale() { return scale_; }
		PolyPath *AddChild(const Path &path) override
		{
			int error_code = 0;
			auto p = std::make_unique<PolyPath>(this);
			auto *result = p.get();
			childs_.emplace_back(std::move(p));

			result->polygon_ = ScalePath(path, scale_, error_code);
			return result;
		}

		void Clear() override
		{
			childs_.resize(0);
		}

		size_t Count() const override
		{
			return childs_.size();
		}

		const Path &Polygon() const { return polygon_; };

		double Area() const
		{
			return std::accumulate(childs_.begin(), childs_.end(), Clipper2Lib::Area(polygon_), [](double a, const auto &child)
								   { return a + child->Area(); });
		}
	};

	class Clipper : public ClipperBase
	{
	private:
		double scale_ = 1.0, invScale_ = 1.0;
		void BuildPaths(Paths &solutionClosed, Paths *solutionOpen);
		void BuildTree(PolyPath &polytree, Paths &open_paths);

	public:
		explicit Clipper(int precision = 7) : ClipperBase()
		{
			CheckPrecision(precision, error_code_);
			// to optimize scaling / descaling precision
			// set the scale to a power of double's radix (2) (#25)
			scale_ = std::pow(std::numeric_limits<double>::radix, std::ilogb(std::pow(10, precision)) + 1);
			invScale_ = 1 / scale_;
		}

		void AddSubject(const Paths &subjects)
		{
			AddPaths(ScalePaths(subjects, scale_, error_code_), PathType::Subject, false);
		}

		void AddOpenSubject(const Paths &open_subjects)
		{
			AddPaths(ScalePaths(open_subjects, scale_, error_code_), PathType::Subject, true);
		}

		void AddClip(const Paths &clips)
		{
			AddPaths(ScalePaths(clips, scale_, error_code_), PathType::Clip, false);
		}

		bool Execute(ClipType clip_type, FillRule fill_rule, Paths &closed_paths)
		{
			Paths dummy;
			return Execute(clip_type, fill_rule, closed_paths, dummy);
		}

		bool Execute(ClipType clip_type,
					 FillRule fill_rule, Paths &closed_paths, Paths &open_paths)
		{
			if (ExecuteInternal(clip_type, fill_rule, false))
			{
				BuildPaths(closed_paths, &open_paths);
			}
			CleanUp();
			return succeeded_;
		}

		bool Execute(ClipType clip_type, FillRule fill_rule, PolyTree &polytree)
		{
			Paths dummy;
			return Execute(clip_type, fill_rule, polytree, dummy);
		}

		bool Execute(ClipType clip_type,
					 FillRule fill_rule, PolyTree &polytree, Paths &open_paths)
		{
			if (ExecuteInternal(clip_type, fill_rule, true))
			{
				polytree.Clear();
				polytree.SetScale(invScale_);
				open_paths.clear();
				BuildTree(polytree, open_paths);
			}
			CleanUp();
			return succeeded_;
		}
	};

} // namespace
