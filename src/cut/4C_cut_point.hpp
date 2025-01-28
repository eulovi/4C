// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_POINT_HPP
#define FOUR_C_CUT_POINT_HPP

#include "4C_config.hpp"

#include "4C_cut_tolerance.hpp"
#include "4C_cut_utils.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_clnwrapper.hpp"

#include <fstream>
#include <map>
#include <memory>

FOUR_C_NAMESPACE_OPEN


namespace Cut
{
  namespace Impl
  {
    class PointLineFilter
    {
     public:
      virtual ~PointLineFilter() = default;
      virtual bool operator()(Line* line) const = 0;
    };

  }  // namespace Impl

  /// Points in a mesh. Both nodal points and cut points.
  /*!
   * A point knows its position with respect to the cut interface.
   */
  class Point
  {
   public:
    enum PointPosition
    {
      undecided = 0,
      oncutsurface = -1,
      inside = -2,
      outside = -3
    };

    /// translate PointPosition enumerator to string
    static inline std::string point_position_to_string(const enum PointPosition& pos)
    {
      switch (pos)
      {
        case undecided:
          return "undecided";
        case oncutsurface:
          return "oncutsurface";
        case inside:
          return "inside";
        case outside:
          return "outside";
        default:
          break;
      }
      return "Unknown PointPosition enum";
    };

    static Point* new_point(Mesh& mesh, const double* x, double t, Edge* cut_edge, Side* cut_side,
        double tolerance = 0.0);

    static Point* insert_cut(Edge* cut_edge, Side* cut_side, Node* n);

    /// constructor
    Point(unsigned pid, double tolerance);

    /// destructor
    virtual ~Point() = default;

    int id() const { return pid_; }

    const double* x() const { return x_; }

    PointPosition position() const { return position_; }

    void position(PointPosition p);

    /*! \brief Add this edge to the list of edges that are cut by this point */
    void add_edge(Edge* e);

    /*! \brief Add this side to the list of sides that are cut by this point */
    void add_side(Side* s);

    /*! \brief Add this element to the list of elements that are cut by this point */
    void add_element(Element* e)
    {
      if (e != nullptr)
      {
        cut_elements_.insert(e);
      }
    }

    /// removing element from the list of elements cut by this point (can be used when position
    /// was determined wrongly)
    void remove_element(Element* e) { cut_elements_.erase(e); }

    /// adding pair of side-edge to container of intersected, does not insert if such pair already
    /// exists
    void add_pair(Side* side, Edge* edge, const std::pair<Side*, Edge*>& original_pair);

    /// default call
    void add_pair(Side* side, Edge* edge);

    /// called during merging
    void add_pair(const std::pair<Side*, Edge*>& pair, const std::pair<Side*, Edge*>& original_pair,
        Point* source);

    /// Add all topological connections of intersection of this edge and other edge ( all
    /// necessery pairs, etc)
    void add_edge_intersection(Edge* first, Edge* second,
        const std::pair<Side*, Edge*>& original_cut_pair, const std::string& extra_msg = "");

    void add_edge_intersection(Edge* first, Edge* second, Side* original_side, Edge* original_edge,
        const std::string& extra_msg = "");

    /// returns true if this point is created by this pair of side and edge intersection
    bool is_cut(Side* side, Edge* edge);

    /*! \brief Returns true if the edge is cut by this point */
    bool is_cut(Edge* e) { return cut_edges_.count(e) > 0; }

    /*! \brief Returns true if the facet is cut by this point */
    bool is_cut(Facet* f) { return facets_.count(f) > 0; }

    /*! \brief Returns true if the side is cut by this point */
    bool is_cut(Side* s) { return cut_sides_.count(s) > 0; }

    /*! \brief Returns true if this point cuts the element */
    bool is_cut(Element* s) { return cut_elements_.count(s) > 0; }

    /// Checks is this cut point is the result of s1 x s2 by checking cut_pairs
    bool is_cut(Side* s1, Side* s2);

    /*
     * Return if this point associated with a facet which is either a cut-side or a marked side
     *  Needed for the tessellation.
     *  TODO: Look into this
     *  winter 02/2017
     */
    bool has_associated_boundary_cell_facet();

    /*! \brief Get the coordinates of this point */
    virtual void coordinates(double* x) const = 0;

    void register_entity(Line* line) { lines_.insert(line); }

    void register_entity(Facet* facet) { facets_.insert(facet); }

    /*!
    \brief Identifies the edges that are cut by considered point and given point
    */
    void common_edge(Point* other, plain_edge_set& edges);

    /*!
    \brief Identifies the sides that are cut by considered point and given point
    */
    void common_side(Point* other, plain_side_set& sides);

    Line* common_line(Point* other);

    Line* cut_line(Side* side, bool unique = true);

    Line* cut_line(Line* line, Side* side, Element* element);

    Line* cut_line(const Impl::PointLineFilter& filter, bool unique = true);

    Line* cut_line(Line* line, const Impl::PointLineFilter& filter, bool unique = true);

    void cut_lines(const Impl::PointLineFilter& filter, plain_line_set& cut_lines);

    void cut_lines(Side* side, plain_line_set& cut_lines);

    // std::vector<Edge*> CutEdges( Point * other );

    const plain_edge_set& cut_edges() { return cut_edges_; }

    const plain_side_set& cut_sides() { return cut_sides_; }

    const std::set<std::pair<Side*, Edge*>>& cut_pairs() { return cut_pairs_; }

    Side* cut_side(Side* side, Point* other);

    void cut_edge(Side* side, Line* other_line, std::vector<Edge*>& matches);

    void print(std::ostream& stream = std::cout) const
    {
      stream << "(" << pid_ << "; " << std::setprecision(16) << x_[0] << ","
             << std::setprecision(16) << x_[1] << "," << std::setprecision(16) << x_[2] << ")";
    }

    void plot(std::ostream& f, int nid = -50) const
    {
      f << std::setprecision(25) << x_[0] << " " << std::setprecision(25) << x_[1] << " "
        << std::setprecision(25) << x_[2] << " # " << pid_ << " " << nid;
      dump_doubles(f, x(), 3);
      f << "\n";
    }

    /** \brief check if the given point lies on the edge and return the
     *  corresponding local coordinate of this point on the edge
     *
     *  Return the parameter space coordinate \f$t\f$ of this point \f$x\f$ on
     *  the given \c edge with begin node \f$\bar{x}_{0}\f$ and the end point
     *  \f$\bar{x}_{1}\f$:
     *
     *  x = \sum\limits_{i=0}^{1} N_{i}(t) \bar{x}_{i}
     *
     *  \param edge (in) : edge which we want to check */
    double t(Edge* edge);

    /// Same as above.
    /// But it can happen that origanl intersection point is merged into different location.
    /// However, we want to obtain position of the cut_point
    /// on the edge based on the real intersection( before mergin),to be consistest with rest of
    /// the code
    double t(Edge* edge, const Core::LinAlg::Matrix<3, 1>& coord);

    // void t( Edge* edge, double pos ) { t_[edge] = pos; }

    unsigned pid() const { return pid_; }

    bool nodal_point(const std::vector<Node*>& nodes) const;

    // checking if the cut_point is close enough to the nodal point
    bool almost_nodal_point(const std::vector<Node*>& nodes, double tolerance = 1e-10) const;

    Node* cut_node();

    void intersection(plain_edge_set& edges);

    void intersection(plain_side_set& sides);

    void intersection(plain_facet_set& facets);

    void intersection(plain_element_set& elements);

    const plain_element_set& elements() const { return cut_elements_; }

    const plain_facet_set& facets() const { return facets_; }

    Edge* common_cut_edge(Side* side);

    /// erase all cut pairs containing this side
    void erased_containing_cut_pairs(Side* cutside);

    /// erased all cut pairs containing this edge
    void erased_containing_cut_pairs(Edge* cutsideedge);
    /// Erase the cutside from this point because it is deleted in the selfcut
    void erase_cut_side(Side* cutside) { cut_sides_.erase(cutside); }

    /// Erase the cutsideedge from this point because it is deleted in the selfcut
    void erase_cut_side_edge(Edge* cutsideedge) { cut_edges_.erase(cutsideedge); }

    /// Return actual point tolerance
    double tolerance() { return tol_; }

    /// Set new Tolerance (in case points are merged this is necessary)
    void enlarge_tolerance(double newtol)
    {
      if (newtol > tol_)
      {
        tol_ = newtol;
      }
    }

    /// Removes information that this point was created by cut from everywhere (e.g. sides,
    /// edges, element)
    void remove_connectivity_info();

    /// Carefully replace this points by another at the stage before/during creation of cut_lines
    void replace(Point* p);

    void remove_edge(Edge* edge);

    void remove_side(Side* side);

    void remove_pair(std::pair<Side*, Edge*>& pair) { cut_pairs_.erase(pair); };

    // move point slighlty (used for merging)
    virtual void move_point(const double* new_coord) = 0;

    // Merge point into anothe r
    void merge(Point* dest);

    /// Dump all information of how this point was created, edges x sides intersections. Useful to
    /// debug reasons crashing
    void dump_connectivity_info();

    /// Get the Merged Points
    std::vector<Point*> get_merged_points() { return merged_points_; };

#if CUT_CREATION_INFO
    // Add merged pair. If the coord = nullptr, it means it was inserted without actual
    // computation of the coordinate, solely on topology e.g. in the insert_cut
    void AddMergedPair(const std::pair<Side*, Edge*>& inter, Core::LinAlg::Matrix<3, 1>* coord);

    const std::map<std::pair<Side*, Edge*>, std::pair<Core::LinAlg::Matrix<3, 1>, bool>>&
    GetMergedPointsInfo()
    {
      return merged_points_info_;
    };

    // Add info about based on which side and edge this point was created (with optional message)
    void AddCreationInfo(const std::pair<Side*, Edge*>& cut_pair, const std::string& info);

    // Add creation info, that does not depended on particular cut pair
    void AddAdditionalCreationInfo(const std::string& info);

    const std::string& GetAdditionalCreationInfo();

    // returns empty string if there is no info for this pair
    const std::string& GetCreationInfo(const std::pair<Side*, Edge*>& pair)
    {
      return creation_info_[pair];
    };


    // for particular cut pair gives info if it is really cut or added later on during topology
    // identification (e.g. because of neigboring edges )
    bool IsReallyCut(const std::pair<Side*, Edge*>& pair)
    {
      std::map<std::pair<Side*, Edge*>, std::pair<std::pair<Side*, Edge*>, Point*>>::iterator it =
          cut_pairs_info_.find(pair);
      if (it != cut_pairs_info_.end())
      {
        std::pair<Side*, Edge*> found_pair = it->second.first;
        return (found_pair == it->first);
      }
      else
      {
        FOUR_C_THROW("Following cut_pair does not exist!");
      }
      return false;
    }

    const std::pair<Side*, Edge*>& AddedFrom(const std::pair<Side*, Edge*>& resultant_pair)
    {
      std::map<std::pair<Side*, Edge*>, std::pair<std::pair<Side*, Edge*>, Point*>>::iterator it =
          cut_pairs_info_.find(resultant_pair);
      if (it != cut_pairs_info_.end())
      {
        return (it->second.first);
      }
      else
      {
        FOUR_C_THROW("Following cut pair does not exists!");
      }
      // should not reach here, just to make compiler happy
      return cut_pairs_info_.begin()->second.first;
    }

#endif

   protected:
    // internal stored point coordinates (always of dimension 3!)
    double x_[3];

   private:
    // find the unique line that matches the filter criteria
    template <class Filter>
    Line* find(Filter& filter, bool unique = true);

    template <class Filter>
    void find(Filter& filter, plain_line_set& cut_lines);

    unsigned pid_;
    PointPosition position_;

    std::map<Edge*, double> t_;

    // precicion of the point
    double tol_;

    plain_edge_set cut_edges_;
    plain_side_set cut_sides_;  // a set with all sides
    plain_element_set cut_elements_;

    plain_line_set lines_;
    plain_facet_set facets_;

    // stores cut pairs of side and edge intersections for this point in associative containers
    // for easy lookup and insert
    std::set<std::pair<Side*, Edge*>> cut_pairs_;

    /// points that were merged here during replacement
    std::vector<Point*> merged_points_;

    /// is this point currently in a process of position definition?
    bool is_point_pos_processed_ = false;

#if CUT_CREATION_INFO

    // Map that indicates whether cut pairs whether originally (during point creation) or later
    // (during neighbors identification), key is the pair that was added, value is the original
    // pair that result in this add. if they are same it was added directly in intersection, if
    // different - during topology identification Point represent origin source point, that was
    // merged into this one
    std::map<std::pair<Side*, Edge*>, std::pair<std::pair<Side*, Edge*>, Point*>> cut_pairs_info_;

    // Map that indicates what intersections (keys of the map)  merged coordinates in the the
    // first part of value pair to the Point* of second part of value parit.Indicates what
    // coordinates were merged in this particular point, due to closeness
    std::map<std::pair<Side*, Edge*>, std::pair<Core::LinAlg::Matrix<3, 1>, bool>>
        merged_points_info_;

    // Displays why particular cut pair is intersected in this point
    std::map<std::pair<Side*, Edge*>, std::string> creation_info_;

    // additional creation info that might be cut pair independent
    std::string additional_creation_info_;

    // Indicates actuall real points, that were merged here (due to mergin in pointpool)
    std::vector<Point*> real_merged_points_;

    // point that this point got merged into
    Point* merged_to_;

#endif
  };  // class Point

  /*--------------------------------------------------------------------------*/
  /* \class ConcretePoint
   *
   * Necessary due to the different problem dimensions. Actually only the
   * constructor and the coordinates() function change, since the stored
   * point coordinate stays 3. If probDim is smaller than three, the last
   * coordinates are set to zero.
   *
   * \author hiermeier
   * \date 11/16 */
  template <unsigned prob_dim>
  class ConcretePoint : public Point
  {
   public:
    /// constructor
    ConcretePoint(unsigned pid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance)
        : Point(pid, tolerance)
    {
      // copy the given coordinates into the base class member variable
      std::copy(x, x + prob_dim, this->x_);

      // set the remaining entries to zero
      std::fill(this->x_ + prob_dim, this->x_ + 3, 0.0);

      // The x_ coordinate must be set, before these methods are called!
      if (cut_edge != nullptr)
      {
        add_edge(cut_edge);
      }
      if (cut_side != nullptr)
      {
        add_side(cut_side);
      }

      if ((cut_side != nullptr) && (cut_edge != nullptr)) this->add_pair(cut_side, cut_edge);
    };

    /** \brief Get the coordinates of this point */
    void coordinates(double* x) const override { std::copy(this->x_, this->x_ + prob_dim, x); };

    void move_point(const double* new_coord) override;

  };  // class ConcretePoint

  /*--------------------------------------------------------------------------*/
  /** \class PointFactory
   *
   *  This class creates a ConcretePoint object with correct problem dimension.
   *  Don't call this class directly! Use the non-member function instead.
   *
   *  \date 11/16
   *  \author hiermeier */
  class PointFactory
  {
   public:
    PointFactory(){};

    // non-member function to create a concrete point of desired dimension
    std::shared_ptr<Cut::Point> create_point(unsigned pid, const double* x, Edge* cut_edge,
        Side* cut_side, double tolerance, int probdim) const
    {
      std::shared_ptr<Cut::Point> point = nullptr;

      switch (probdim)
      {
        case 2:
          point = std::make_shared<ConcretePoint<2>>(pid, x, cut_edge, cut_side, tolerance);
          break;
        case 3:
          point = std::make_shared<ConcretePoint<3>>(pid, x, cut_edge, cut_side, tolerance);
          break;
        default:
          FOUR_C_THROW("Unsupported problem dimension! (probdim=%d)", probdim);
          break;
      }
      return point;
    };
  };  // class PointFactory

  // non-member function to create a concrete point of desired dimension
  std::shared_ptr<Cut::Point> create_point(
      unsigned pid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance);

  inline int entity_id(const Point& p) { return p.pid(); }

  /// id based comparison for sorting and searching
  template <class T>
  class EntityIdLess
  {
   public:
    bool operator()(const T& p1, const T& p2) const { return EntityId(p1) < EntityId(p2); }

    bool operator()(const T* p1, const T* p2) const { return entity_id(*p1) < entity_id(*p2); }

    bool operator()(const std::shared_ptr<T> p1, const std::shared_ptr<T> p2) const
    {
      (void)p1;
      (void)p2;
      return entity_id(*p1) < entity_id(*p2);
    }
  };

  typedef EntityIdLess<Point> PointPidLess;

  /// position on edge based comparison for sorting and searching
  class PointPositionLess
  {
   public:
    explicit PointPositionLess(Edge* edge) : edge_(edge) {}

    bool operator()(Point& p1, Point& p2) const
    {
      if (not p1.is_cut(edge_) or not p2.is_cut(edge_))
        FOUR_C_THROW("point position compare only on cut edges");
      return p1.t(edge_) < p2.t(edge_);
    }

    bool operator()(Point* p1, Point* p2) const
    {
      if (not p1->is_cut(edge_) or not p2->is_cut(edge_))
        FOUR_C_THROW("point position compare only on cut edges");
      return p1->t(edge_) < p2->t(edge_);
    }

   private:
    Edge* edge_;
  };

#ifdef CUT_USE_SORTED_VECTOR
  typedef SortedVector<Point*, true, PointPidLess> PointSet;
  typedef SortedVector<Point*, true, PointPositionLess> PointPositionSet;

  typedef SortedVector<std::shared_ptr<Point>, true, PointPidLess> RCPPointSet;

#else
  typedef std::set<Point*, PointPidLess> PointSet;
  typedef std::set<Point*, PointPositionLess> PointPositionSet;

  typedef std::set<std::shared_ptr<Point>, PointPidLess> RCPPointSet;

#endif

  namespace Impl
  {
    class ExcludeLineFilter : public PointLineFilter
    {
     public:
      ExcludeLineFilter(Line* line, const PointLineFilter& filter) : line_(line), filter_(filter) {}

      bool operator()(Line* line) const override { return line != line_ and filter_(line); }

     private:
      Line* line_;
      const PointLineFilter& filter_;
    };

    class SideCutFilter : public PointLineFilter
    {
     public:
      SideCutFilter(Side* side) : side_(side) {}

      bool operator()(Line* line) const override;

     private:
      Side* side_;
    };

    class SideElementCutFilter : public PointLineFilter
    {
     public:
      SideElementCutFilter(Side* side, Element* element) : side_(side), element_(element) {}

      bool operator()(Line* line) const override;

     private:
      Side* side_;
      Element* element_;
    };

    class SideSideCutFilter : public PointLineFilter
    {
     public:
      SideSideCutFilter(Side* side1, Side* side2) : side1_(side1), side2_(side2) {}

      bool operator()(Line* line) const override;

     private:
      Side* side1_;
      Side* side2_;
    };

  }  // namespace Impl

  /** \brief Find if the points in element share a common element,
   *  (i.e. do all points lie on the same element?)
   *
   *  A point knows which cut_elements it is associated to.
   *
   *  \author hiermeier \date 11/16 */
  void find_common_elements(const std::vector<Point*>& element, plain_element_set& elements);

  /// Find if the points in side (which is actually a tet...) share a common side,
  ///   (i.e. do all points lie on the same surface?)
  /// A point knows which cut_sides it is associated to.
  inline void find_common_sides(const std::vector<Point*>& side, plain_side_set& sides)
  {
    std::vector<Point*>::const_iterator is = side.begin();
    // Get the sides this point cuts
    // A cut-side which is a LevelSet-side is empty but is still present.
    sides = (*is)->cut_sides();
    for (++is; is != side.end(); ++is)
    {
      Point* p = *is;
      p->intersection(sides);
      if (sides.size() == 0)
      {
        break;
      }
    }
  }

  inline void find_common_sides(Point* p1, Point* p2, Point* p3, plain_side_set& sides)
  {
    sides = p1->cut_sides();
    p2->intersection(sides);
    p3->intersection(sides);
  }

  inline void find_common_sides(Point* p1, Point* p2, Point* p3, Point* p4, plain_side_set& sides)
  {
    sides = p1->cut_sides();
    p2->intersection(sides);
    p3->intersection(sides);
    p4->intersection(sides);
  }

  /// Find distance between points
  template <unsigned int prob_dim>
  double distance_between_points(const Core::LinAlg::Matrix<prob_dim, 1>& coord_a,
      const Core::LinAlg::Matrix<prob_dim, 1>& coord_b)
  {
    Core::LinAlg::Matrix<prob_dim, 1> diff;
    diff.update(1, coord_a, -1, coord_b);
    return diff.norm2();
  }

  /// Find distance between points
  double distance_between_points(Point* p1, Point* p2);

  /// Find distance between points
  double distance_between_points(Point* p1, const Core::LinAlg::Matrix<3, 1>& coord_b);

  bool is_cut_position_unchanged(Point::PointPosition position, Point::PointPosition pos);

}  // namespace Cut



std::ostream& operator<<(std::ostream& stream, Cut::Point& point);

std::ostream& operator<<(std::ostream& stream, Cut::Point* point);

FOUR_C_NAMESPACE_CLOSE

#endif
