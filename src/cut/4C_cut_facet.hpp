// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_FACET_HPP
#define FOUR_C_CUT_FACET_HPP

#include "4C_config.hpp"

#include "4C_cut_line.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_io_pstream.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Cut
{
  class BoundaryCell;
  class VolumeCell;
  class Element;
  class Mesh;
  class Facet;

  /*!
  \brief Class to handle surfaces of arbitrary shape, defined by its points corner points
   */
  class Facet
  {
   public:
    Facet(Mesh& mesh, const std::vector<Point*>& points, Side* side, bool cutsurface);

    void register_entity(VolumeCell* cell);

    void disconnect_volume(VolumeCell* cell);

    inline void print() const { print(std::cout); }

    void print(std::ostream& stream) const;

    /// Print only point IDs of a facet.
    void print_point_ids()
    {
      for (std::vector<Point*>::iterator i = points_.begin(); i != points_.end(); ++i)
      {
        Point& p = **i;
        std::cout << p.pid() << " ";
      }
      std::cout << "\n";
    }



    /** TRUE, if the \c parentside_ has a ID greater than -1 and is thus no
     *  \c element_side, i.e. side of the background mesh. */
    bool on_cut_side() const;

    /*!
    \brief Return true if this facet is on a marked side from the background mesh
     */
    bool on_marked_background_side() const;

    /* Check if the Facet belongs to a side which is either cut OR marked
     * -> i.e. it should create boundary cells.
     *  */
    bool on_boundary_cell_side() const;

    /*!
    \brief Returns the parent side Id from which the facet is created
     */
    int side_id() const;

    Side* parent_side() const { return parentside_; }

    void coordinates(double* x);

    void corner_coordinates(double* x);

    void get_all_points(Mesh& mesh, PointSet& cut_points, bool dotriangulate = false);

    void add_hole(Facet* hole);

    /** \brief set the given side as parentside and set the position as well */
    void exchange_side(Side* side, bool cutsurface)
    {
      parentside_ = side;
      if (cutsurface)
      {
        position(Point::oncutsurface);
        for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
        {
          Point* p = *i;
          p->position(Point::oncutsurface);
        }
      }
    }

    bool equals(const std::vector<Point*>& facet_points) { return equals(points_, facet_points); }

    bool equals(Core::FE::CellType distype);

    bool corner_equals(const std::vector<Point*>& facet_points)
    {
      return equals(corner_points_, facet_points);
    }

    /*! \brief Check whether the parent side is a cut side */
    bool is_cut_side(Side* side);

    Point::PointPosition position() const { return position_; }

    void position(Point::PointPosition p);

    void get_lines(std::map<std::pair<Point*, Point*>, plain_facet_set>& lines);


    void get_lines(const std::vector<Point*>& points,
        std::map<std::pair<Point*, Point*>, plain_facet_set>& lines);

    bool is_line(Point* p1, Point* p2);

    bool contains(Point* p) const;


    /** \brief Check if the given volume cell is equal to one of the already
     *         stored volume cells in this facet.
     *
     *  \author  hiermeier \date 12/16 */
    bool contains(const plain_facet_set& vcell) const;

    bool contains(const std::vector<Point*>& side) const;

    bool contains_some(const std::vector<Point*>& side) const;

    bool touches(Facet* f);

    /*!
    \brief If this Facet has a CommonEdge with another facet, based on this edge the point
    ordering is checked
    */
    bool have_consistant_normal(Facet* f,  // f ... facetpointer to facet to compare with!
        bool& result);  // result == true --> normal points in the same direction!

    VolumeCell* neighbor(VolumeCell* cell);

    void neighbors(Point* p, const plain_volumecell_set& cells, const plain_volumecell_set& done,
        plain_volumecell_set& connected, plain_element_set& elements);

    void neighbors(Point* p, const plain_volumecell_set& cells, const plain_volumecell_set& done,
        plain_volumecell_set& connected);

    const std::vector<Point*>& points() const { return points_; }

    /*!
    \brief Get the corner points of the facet in global coordinates
     */
    const std::vector<Point*>& corner_points() const { return corner_points_; }

    /*!
    \brief Get the corner points of the facet in element local coordinates. Used in Moment fitting
    method
     */
    void corner_points_local(
        Element* elem1, std::vector<std::vector<double>>& cornersLocal, bool shadow = false);

    /*!
    \brief Return the global coordinates all of its corner points in order
     */
    std::vector<std::vector<double>> corner_points_global(Element* elem1, bool shadow = false);

    /*!
    \brief Get the triangulated sides of this facet
     */
    const std::vector<std::vector<Point*>>& triangulation() const { return triangulation_; }

    /*!
    \brief Get all the triangulated points in the specified pointset
     */
    void triangulation_points(PointSet& points);

    void all_points(PointSet& points)
    {
      if (is_triangulated())
      {
        triangulation_points(points);
      }
      else
      {
        std::copy(points_.begin(), points_.end(), std::inserter(points, points.begin()));
      }
    }

    /// Create new point1 boundary cell associated with this facet
    void new_point1_cell(Mesh& mesh, VolumeCell* volume, const std::vector<Point*>& points,
        plain_boundarycell_set& bcells);

    /// Create new point1 boundary cell associated with this facet
    void new_line2_cell(Mesh& mesh, VolumeCell* volume, const std::vector<Point*>& points,
        plain_boundarycell_set& bcells);

    /*!
    \brief Create new tri3 boundary cell associated with this facet
     */
    void new_tri3_cell(Mesh& mesh, VolumeCell* volume, const std::vector<Point*>& points,
        plain_boundarycell_set& bcells);

    /*!
    \brief Create new quad4 boundary cell associated with this facet
     */
    void new_quad4_cell(Mesh& mesh, VolumeCell* volume, const std::vector<Point*>& points,
        plain_boundarycell_set& bcells);

    /*!
    \brief Create new arbitrary boundary cell associated with this facet. These cells are to be
    dealt with when moment fitting is used for boun.cell integration
     */
    void new_arbitrary_cell(Mesh& mesh, VolumeCell* volume, const std::vector<Point*>& points,
        plain_boundarycell_set& bcells, const Core::FE::GaussIntegration& gp,
        const Core::LinAlg::Matrix<3, 1>& normal);

    /// Get the BoundaryCells created on this facet
    void get_boundary_cells(plain_boundarycell_set& bcells);

    void test_facet_area(double tolerance, bool istetmeshintersection = false);

    bool is_triangle(const std::vector<Point*>& tri) const;

    /*!
    \brief Check whether the facet is already triangulated
     */
    bool is_triangulated() const { return triangulation_.size() > 0; }

    /*!
    \brief Check whether the given vector of points is a triangulation of this facet
     */
    bool is_triangulated_side(const std::vector<Point*>& tri) const;

    bool has_holes() const { return holes_.size() > 0; }

    const plain_facet_set& holes() const { return holes_; }

    unsigned num_points();

    const plain_volumecell_set& cells() const { return cells_; }

    Point* other_point(Point* p1, Point* p2);

    /*!
    \brief Triangulate the facet. This happens implicitly if Tessellation is used. This simply
    triangulates the facet any may not give outward normal for the resulting cells
     */
    void do_triangulation(Mesh& mesh, const std::vector<Point*>& points)
    {
      create_triangulation(mesh, points);
    }

    /*!
    \brief check whether facet is already split
     */
    bool is_facet_split() const { return split_cells_.size() > 0; }

    /*!
    \brief split the facet into a number of tri and quad. Reduced number of Gauss points when
    facet is split instead of triangulated
     */
    void split_facet(const std::vector<Point*>& facetpts);

    /*!
    \brief Get the triangulated sides of this facet
     */
    const std::vector<std::vector<Point*>>& get_split_cells() const { return split_cells_; }

    bool is_planar(Mesh& mesh, const std::vector<Point*>& points);

    /// Do the facets share the same CutSide?
    bool share_same_cut_side(Facet* f);

    /// Return true is the facet is convex shaped
    bool is_convex();

    /// Belongs to a LevelSetSide
    bool belongs_to_level_set_side();

   private:
    Facet(const Facet&);
    Facet& operator=(const Facet&);

    bool is_planar(Mesh& mesh, bool dotriangulate);

    void create_triangulation(Mesh& mesh, const std::vector<Point*>& points);

    void get_nodal_ids(Mesh& mesh, const std::vector<Point*>& points, std::vector<int>& nids);

    unsigned normal(const std::vector<Point*>& points, Core::LinAlg::Matrix<3, 1>& x1,
        Core::LinAlg::Matrix<3, 1>& x2, Core::LinAlg::Matrix<3, 1>& x3,
        Core::LinAlg::Matrix<3, 1>& b1, Core::LinAlg::Matrix<3, 1>& b2,
        Core::LinAlg::Matrix<3, 1>& b3);

    void find_corner_points();

    bool is_line(const std::vector<Point*>& points, Point* p1, Point* p2);

    bool equals(const std::vector<Point*>& my_points, const std::vector<Point*>& facet_points);

    std::vector<Point*> points_;

    std::vector<Point*> corner_points_;

    plain_facet_set holes_;

    std::vector<std::vector<Point*>> triangulation_;

    std::vector<std::vector<Point*>> split_cells_;

    Side* parentside_;

    bool planar_;

    bool planar_known_;

    Point::PointPosition position_;

    plain_volumecell_set cells_;

    //! already cheacked whether the facet is planar?
    bool is_planar_computed_;

    //! whether this facet lie on a plane?
    bool is_planar_;

    /// is this facet currently in a process of position definition?
    bool is_facet_pos_processed_ = false;
  };

  template <class T>
  inline Facet* find_facet(const T& facets, const std::vector<Point*>& side)
  {
    Facet* found = nullptr;
    for (typename T::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f->corner_equals(side))
      {
        if (found == nullptr)
        {
          found = f;
        }
        else
        {
          FOUR_C_THROW("not unique");
        }
      }
    }
    return found;
  }

  // inline int EntityId( const Facet & f ) { return f.InternalId(); }
  // inline int EntityId( const Facet & f ) { return reinterpret_cast<int>( &f ); }

  inline void remove_nonmatching_triangulated_facets(
      const std::vector<Point*>& side, plain_facet_set& facets)
  {
    if (side.size() == 3)
    {
      for (plain_facet_set::iterator i = facets.begin(); i != facets.end();)
      {
        Facet* f = *i;
        if (f->is_triangulated())
        {
          if (not f->is_triangulated_side(side))
          {
            set_erase(facets, i);
          }
          else
          {
            ++i;
          }
        }
        else
        {
          ++i;
        }
      }
    }
  }

  inline void find_common_facets(const std::vector<Point*>& side, plain_facet_set& facets)
  {
    std::vector<Point*>::const_iterator is = side.begin();
    facets = (*is)->facets();
    for (++is; is != side.end(); ++is)
    {
      Point* p = *is;
      p->intersection(facets);
      if (facets.size() == 0)
      {
        break;
      }
    }
    // This is probably an unnecessary call as side here is a tet, i.e. side.size() == 4.
    if (side.size() == 3) FOUR_C_THROW("The TET is degenerate! It does not contain 4 points!");
    // Might be able to remove this call requires side.size()==3
    remove_nonmatching_triangulated_facets(side, facets);
  }

  inline void find_common_facets(Point* p1, Point* p2, Point* p3, plain_facet_set& facets)
  {
    facets = p1->facets();
    p2->intersection(facets);
    p3->intersection(facets);

    std::vector<Point*> side(3);
    side[0] = p1;
    side[1] = p2;
    side[2] = p3;
    remove_nonmatching_triangulated_facets(side, facets);
  }

  inline void find_common_facets(
      Point* p1, Point* p2, Point* p3, Point* p4, plain_facet_set& facets)
  {
    facets = p1->facets();
    p2->intersection(facets);
    p3->intersection(facets);
    p4->intersection(facets);
  }

}  // namespace Cut


std::ostream& operator<<(std::ostream& stream, Cut::Facet& f);

FOUR_C_NAMESPACE_CLOSE

#endif
