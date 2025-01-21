// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_ELEMENTHANDLE_HPP
#define FOUR_C_CUT_ELEMENTHANDLE_HPP


#include "4C_config.hpp"

#include "4C_cut_element.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_utils_pairedvector.hpp"

FOUR_C_NAMESPACE_OPEN



namespace Cut
{
  class Point;
  class Node;
  class Element;
  class Mesh;


  /*!
  \brief Outside world interface to element. Converts quadratic to linear element. This provides
  the Gaussian rules generated from the cut
   */
  //--------------------------------------------------------------------------//
  // the ElementHandle base class
  //--------------------------------------------------------------------------//
  class ElementHandle
  {
   public:
    virtual ~ElementHandle() = default;
    //--------------------------------------------------------------------------//
    //! @name access to cut/intersection status of the element

    /*!
   \brief returns true in case that any cut-side cut with the element produces cut points,
          i.e. also for touched cases (at points, edges or sides),
          or when an element side has more than one facet or is touched by fully/partially
          by the cut side
     */
    virtual bool is_cut() = 0;

    /// return true if the element has more than one volume-cell and therefore is intersected by a
    /// cut-side
    virtual bool is_intersected() = 0;

    //@}

    //--------------------------------------------------------------------------//
    //! @name access to basic element data

    //! get the shape of the element
    virtual Core::FE::CellType shape() = 0;

    //! get the nodes of the element
    virtual const std::vector<Node*>& nodes() = 0;

    //! compute local coordinates of the element for a given point
    const Core::LinAlg::Matrix<3, 1>& local_coordinates(Point* p)
    {
      std::map<Point*, Core::LinAlg::Matrix<3, 1>>::iterator i = local_coordinates_.find(p);
      if (i != local_coordinates_.end())
      {
        return i->second;
      }
      Core::LinAlg::Matrix<3, 1>& rst = local_coordinates_[p];
      Core::LinAlg::Matrix<3, 1> xyz;
      p->coordinates(xyz.data());
      local_coordinates(xyz, rst);
      return rst;
    }

    //! compute local coordinates of the element for given global coordinates
    virtual void local_coordinates(
        const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<3, 1>& rst) = 0;


    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's Cut::Elements (sub-elements)

    virtual void collect_elements(plain_element_set& elements) = 0;

    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's volume-cells

    //! Collect all volume-cells belonging to this elements
    virtual void get_volume_cells(plain_volumecell_set& cells) = 0;

    //! Collect all volume-cells belonging to this element ordered by position
    virtual void collect_volume_cells(
        plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside) = 0;

    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's volume integration specific quantities

    //! get the element's volumetric integration cells (just for Tessellation)
    virtual void get_integration_cells(plain_integrationcell_set& cells) = 0;


    /*!
     \brief  Collect the Gaussian points of all volume-cells belonging to this element in such a
     way that Gaussian rule for every volume-cell can be separated
      */
    void volume_cell_gauss_points(
        plain_volumecell_set& cells, std::vector<Core::FE::GaussIntegration>& intpoints);

    void append_volume_cell_gauss_points_tessellation(
        Core::FE::GaussPointsComposite& gpc, Cut::VolumeCell* vc);

    void append_volume_cell_gauss_points_moment_fitting(
        Core::FE::GaussPointsComposite& gpc, Cut::VolumeCell* vc);

    void append_volume_cell_gauss_points_direct_divergence(
        std::shared_ptr<Core::FE::GaussPointsComposite> gpc, Cut::VolumeCell* vc);


    /*!
     \brief Collect the Gaussian points of all the volume-cells belonging to this element.
            The integration rules over all the volume-cells are connected.
     */
    std::shared_ptr<Core::FE::GaussPointsComposite> gauss_points_connected(
        plain_volumecell_set& cells, VCellGaussPts gausstype);


    //! ...
    virtual int num_volume_cell_gauss_points() = 0;

    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's interface integration specific quantities

    //! get the element's boundary integration cells
    // TODO note which ones???
    virtual void get_boundary_cells(plain_boundarycell_set& bcells) = 0;

    /*!
    \brief Collect the Gauss points of all the boundary-cells belong to this element.
    This is the method used now in the new implementation
     */
    void boundary_cell_gauss_points_lin(
        const std::map<int, std::vector<Cut::BoundaryCell*>>& bcells,
        std::map<int, std::vector<Core::FE::GaussIntegration>>& intpoints,
        const int bc_cubaturedegree);

    //@}


    //--------------------------------------------------------------------------//
    //! @name access to the element's sets of volume-cells, nds-vectors and integration points

    //! get all the element' sets of volume-cells, nds-vectors and integration points
    virtual bool get_cell_sets_dof_sets_gauss_points(std::vector<plain_volumecell_set>& cell_sets,
        std::vector<std::vector<int>>& nds_sets,
        std::vector<std::vector<Core::FE::GaussIntegration>>& intpoints_sets, bool include_inner);

    /*!
    \brief Get the Gauss rule for integration cells.
     *
     * The volume cells outside the integration domain are ignored. Sometimes, after the
     * cut, cells marked as inside cells are actually outside the integration domain. If
     * the parameter integrate_inside_cells == true, the cells marked as inside cells are
     * integrated. Otherwise, the cells marked as outside cells are integrated.
     *
     * @param[out] intpoints_vec                 vector with gauss integration rule
     * @param[in]  integrate_inside_volumecells  boolean value to check if the inside cells should
     * be integrated
     */
    virtual bool get_gauss_rule_integration_cells(
        std::vector<std::vector<Core::FE::GaussIntegration>>& intpoints_sets,
        bool integrate_inside_volumecells);

    //! get the element's sets of volume-cells ordered by inside/outside position
    virtual void volume_cell_sets() = 0;

    //! get all the element' sets of volume-cells and nds-vectors
    void get_volume_cells_dof_sets(std::vector<plain_volumecell_set>& cellsets,
        std::vector<std::vector<int>>& nds_sets, bool include_inner);

    void get_boundary_cell_sets(const std::vector<Cut::Point::PointPosition>& desired_positions,
        std::vector<plain_boundarycell_set>& bcellsets);

    void get_boundary_cell_sets(
        Cut::Point::PointPosition desired_position, std::vector<plain_boundarycell_set>& bcellsets)
    {
      const std::vector<Cut::Point::PointPosition> desired_positions(1, desired_position);
      get_boundary_cell_sets(desired_positions, bcellsets);
    }

   protected:
    /// fill the derived class internal boundary cell sets
    virtual void boundary_cell_set(Point::PointPosition position) = 0;

    /// access the class internal boundary cells sets
    virtual const std::vector<plain_boundarycell_set>& get_boundary_cell_set(
        Point::PointPosition position) = 0;

   public:
    //! get the element's sets of volume-cells with inside position
    virtual const std::vector<plain_volumecell_set>& get_vc_sets_inside() = 0;

    //! get the element's sets of volume-cells with outside position
    virtual const std::vector<plain_volumecell_set>& get_vc_sets_outside() = 0;

    //! get the element's sets of nodal dofset vectors (nds-vectors) with inside position
    std::vector<std::vector<int>>& get_nodal_dof_set_vc_sets_inside()
    {
      return nodaldofset_vc_sets_inside_;
    };

    //! get the element's sets of nodal dofset vectors (nds-vectors) with outside position
    std::vector<std::vector<int>>& get_nodal_dof_set_vc_sets_outside()
    {
      return nodaldofset_vc_sets_outside_;
    };

    //! ...
    std::vector<std::map<int, int>>& get_node_dofset_map_vc_sets_inside_for_communication()
    {
      return vcsets_nid_dofsetnumber_toComm_inside_;
    };

    //! ...
    std::vector<std::map<int, int>>& get_node_dofset_map_vc_sets_outside_for_communication()
    {
      return vcsets_nid_dofsetnumber_toComm_outside_;
    };

    //@}


   private:
    /*!
     \brief Project the integration rule available in the local coordinates of the
     integation-cells to the local coordinates of background element
     */
    template <Core::FE::CellType distype>
    std::shared_ptr<Core::FE::GaussPoints> create_projected(
        const std::vector<Cut::Point*>& cpoints, std::shared_ptr<Core::FE::GaussPoints> gp_ic);

    std::map<Point*, Core::LinAlg::Matrix<3, 1>> local_coordinates_;

   protected:
    /// dof set number of all element nodes, contains the dofset numbers for all nodes of the
    /// superior element (i.e. 20 for hex20 superior element)
    std::vector<std::vector<int>> nodaldofset_vc_sets_inside_;
    std::vector<std::vector<int>> nodaldofset_vc_sets_outside_;

    /// for each set of volumecells a map containing nids and dofsetnumbers that has to be
    /// communicated for the superior element (i.e. 20 for hex20 superior element)
    std::vector<std::map<int, int>> vcsets_nid_dofsetnumber_toComm_inside_;
    std::vector<std::map<int, int>> vcsets_nid_dofsetnumber_toComm_outside_;
  };



  //--------------------------------------------------------------------------//
  // the linear ElementHandle for linear elements
  //--------------------------------------------------------------------------//
  /// linear element handle
  class LinearElementHandle : public ElementHandle
  {
   public:
    LinearElementHandle() : element_(nullptr), cells_set_(false), bcell_sets_(0) {}

    explicit LinearElementHandle(Element* e) : element_(e), cells_set_(false), bcell_sets_(0)
    {
      // set also the parent element Id which is trivial the same as the element Id
      element_->parent_id(e->id());
    }

    //--------------------------------------------------------------------------//
    //! @name access to cut/intersection status of the element

    /*!
     \brief returns true in case that any cut-side cut with the element produces cut points,
            i.e. also for touched cases (at points, edges or sides),
            or when an element side has more than one facet or is touched by fully/partially by
     the cut side
     */
    bool is_cut() override { return element_->is_cut(); }

    /// return true if the element has more than one volume-cell and therefore is intersected by a
    /// cut-side
    bool is_intersected() override { return element_->is_intersected(); }

    //@}


    //--------------------------------------------------------------------------//
    //! @name access to basic element data

    //! get the shape of the element
    Core::FE::CellType shape() override { return element_->shape(); }

    //! get the nodes of the element
    const std::vector<Node*>& nodes() override { return element_->nodes(); }

    //! compute local coordinates of the element for given global coordinates
    void local_coordinates(
        const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<3, 1>& rst) override
    {
      element_->local_coordinates(xyz, rst);
    }

    //@}


    //--------------------------------------------------------------------------//
    //! @name access to the element's Cut::Elements (subelements)

    //! collect all sub-elements
    void collect_elements(plain_element_set& elements) override { elements.insert(element_); }
    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's volume-cells

    //! Collect all volume-cells belonging to this elements
    void get_volume_cells(plain_volumecell_set& cells) override;

    //! Collect all volume-cells belonging to this element ordered by position
    void collect_volume_cells(
        plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside) override;

    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's volume integration specific quantities

    //! get the element's volumetric integration cells (just for Tessellation)
    void get_integration_cells(plain_integrationcell_set& cells) override
    {
      element_->get_integration_cells(cells);
    }

    //! ...
    int num_volume_cell_gauss_points() override { return element_->num_gauss_points(shape()); }

    //@}



    //--------------------------------------------------------------------------//
    //! @name access to the element's interface integration specific quantities

    //! get the element's boundary integration cells
    // TODO note which ones???
    void get_boundary_cells(plain_boundarycell_set& bcells) override
    {
      FOUR_C_THROW(
          "Deprecated function! Use the get_boundary_cell_set( Point::PointPosition ), "
          "instead!");
      element_->get_boundary_cells(bcells);
    }

   protected:
    /** \brief Get the boundary cell sets with given position
     *
     *  \param position (in) : desired position of the boundary cell
     *
     *  \author hiermeier \date 01/17 */
    const std::vector<plain_boundarycell_set>& get_boundary_cell_set(
        Point::PointPosition position) override
    {
      boundary_cell_set(position);
      return bcell_sets_.at(position);
    }

    /** \brief Fill the class internal boundary sets for the given position
     *
     *  This function does a direct return, if the boundary cells have been already
     *  extracted, otherwise we extract the boundary cells from the volume cells
     *  and add them.
     *
     *  \param position (in) : desired position of the boundary cell
     *
     *  \author hiermeier \date 01/17 */
    void boundary_cell_set(Point::PointPosition position) override;

    //@}

   public:
    //--------------------------------------------------------------------------//
    //! @name access to the element's sets of volume-cells, nds-vectors and integration points

    //! get the element's sets of volume-cells ordered by inside/outside position
    void volume_cell_sets() override;

    //! get the element's sets of volume-cells with inside position
    const std::vector<plain_volumecell_set>& get_vc_sets_inside() override
    {
      if (!cells_set_) volume_cell_sets();  // build the volume cell sets
      return vc_sets_inside_;
    };

    //! get the element's sets of volume-cells with outside position
    const std::vector<plain_volumecell_set>& get_vc_sets_outside() override
    {
      if (!cells_set_) volume_cell_sets();  // build the volume cell sets
      return vc_sets_outside_;
    };

    //@}


   private:
    Element* element_;  ///< the unique element in the cut library


   protected:
    bool cells_set_;  ///< have the cells already set?
    std::vector<plain_volumecell_set>
        vc_sets_inside_;  ///< connected sets of volume-cells with inside position
    std::vector<plain_volumecell_set>
        vc_sets_outside_;  ///< connected sets of volume-cells with outside position


    Core::Gen::Pairedvector<Point::PointPosition, std::vector<plain_boundarycell_set>> bcell_sets_;
  };



  //--------------------------------------------------------------------------//
  // the quadratic ElementHandle base class for quadratic elements
  //--------------------------------------------------------------------------//
  /// quadratic element handle
  class QuadraticElementHandle : public ElementHandle
  {
   public:
    QuadraticElementHandle() : cells_connected_(false), connected_bcell_sets_(0) {}

    //--------------------------------------------------------------------------//
    //! @name access to cut/intersection status of the element

    /*!
      \brief returns true in case that any cut-side cut with the element produces cut points,
             i.e. also for touched cases (at points, edges or sides),
             or when an element side has more than one facet or is touched by fully/partially by
      the cut side
     */
    bool is_cut() override;

    /// return true if the element has more than one volume-cell and therefore is intersected by a
    /// cut-side
    bool is_intersected() override;

    //@}


    //--------------------------------------------------------------------------//
    //! @name access to basic element data

    //! get the nodes of the element
    const std::vector<Node*>& nodes() override { return nodes_; }

    //@}


    //--------------------------------------------------------------------------//
    //! @name access to the element's Cut::Elements (subelements)

    //! collect all sub-elements
    void collect_elements(plain_element_set& elements) override
    {
      std::copy(
          subelements_.begin(), subelements_.end(), std::inserter(elements, elements.begin()));
    }
    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's volume-cells

    //! Collect all volume-cells belonging to this elements
    void get_volume_cells(plain_volumecell_set& cells) override;

    //! Collect all volume-cells belonging to this element ordered by position
    void collect_volume_cells(
        plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside) override;

    //! Collect volume-cells belonging to this element with the given position
    void collect_volume_cells(
        Point::PointPosition position, plain_volumecell_set& evolcells_position) const;


    //@}

    //--------------------------------------------------------------------------//
    //! @name access to the element's volume integration specific quantities

    //! get the quadratic element's volumetric integration cells (just for Tessellation)
    void get_integration_cells(plain_integrationcell_set& cells) override;

    //! ...
    int num_volume_cell_gauss_points() override
    {
      int numgp = 0;
      for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
      {
        Element* e = *i;
        numgp += e->num_gauss_points(shape());
      }
      return numgp;
    }

    //@}



    //--------------------------------------------------------------------------//
    //! @name access to the element's interface integration specific quantities

    /// @{

    //! get all the quadratic element's boundary integration cells
    // TODO note which ones???
    void get_boundary_cells(plain_boundarycell_set& bcells) override;

   protected:
    /** \brief Fill the class internal connected boundary sets for the given position
     *
     *  This function does a direct return, if the connected boundary cells have been
     *  already created, otherwise we extract the boundary cells from the connected
     *  volume cells and add them.
     *
     *  \param position (in) : desired position of the boundary cell
     *
     *  \author hiermeier \date 01/17 */
    void boundary_cell_set(Point::PointPosition position) override;

    /** \brief Get the connected boundary cell sets with given position
     *  ( read access )
     *
     *  \param position (in) : desired position of the boundary cell
     *
     *  \author hiermeier \date 01/17 */
    const std::vector<plain_boundarycell_set>& get_boundary_cell_set(
        Point::PointPosition position) override
    {
      boundary_cell_set(position);
      return connected_bcell_sets_.at(position);
    }

   private:
    /** \brief Connect volume cells first and start to create the corresponding
     *  boundary cells afterwards.
     *
     *  \param position (in) : desired position of the boundary cell
     *
     *  \author hiermeier \date 01/17 */
    void connect_boundary_cells(Point::PointPosition position);

    /** \brief After the volume cells have been connected, the connected boundary
     *  cells are build
     *
     *  \param connected_vcell_set (in)  : connected volume cell sets
     *  \param connected_bcell_set (out) : connected boundary cell sets
     *
     *  \author hiermeier \date 01/17 */
    void build_boundary_cell_sets(const std::vector<plain_volumecell_set>& connected_vcell_set,
        std::vector<plain_boundarycell_set>& connected_bcell_set) const;

    /// @}

   public:
    //--------------------------------------------------------------------------//
    //! @name access to the element's sets of volume-cells, nds-vectors and integration points
    /// @{
    //! get the element's sets of volume-cells ordered by inside/outside position
    void volume_cell_sets() override;

    //! get connections/sets of volume-cells between sub-elements ordered by inside position
    const std::vector<plain_volumecell_set>& get_vc_sets_inside() override
    {
      if (!cells_connected_) volume_cell_sets();
      return connected_vc_sets_inside_;
    };

    //! get connections/sets of volume-cells between sub-elements ordered by outside position
    const std::vector<plain_volumecell_set>& get_vc_sets_outside() override
    {
      if (!cells_connected_) volume_cell_sets();
      return connected_vc_sets_outside_;
    };

    /// @}


    //--------------------------------------------------------------------------//
    //! @name specific routines to connect volume-cells between sub-elements which belong to one
    //! connection

    //! connect volume-cells to sets of volume-cells
    virtual void connect_volume_cells();

    //! build sets
    void build_cell_sets(
        plain_volumecell_set& cells_to_connect, std::vector<plain_volumecell_set>& connected_sets);

    //@}


   protected:
    bool cells_connected_;  ///< have the volume-cells been already connected between sub-elements?
    std::vector<plain_volumecell_set>
        connected_vc_sets_inside_;  ///< connected volume-cells with inside position
    std::vector<plain_volumecell_set>
        connected_vc_sets_outside_;  ///< connected volume-cells with outside position

    Core::Gen::Pairedvector<Point::PointPosition, std::vector<plain_boundarycell_set>>
        connected_bcell_sets_;

    std::vector<Element*> subelements_;  ///< the quadratic element's linear sub-elements
    std::vector<Node*> nodes_;           ///< the quadratic element's nodes
  };

  /// hex20 element handle
  class Hex20ElementHandle : public QuadraticElementHandle
  {
   public:
    //! constructor
    Hex20ElementHandle(Mesh& mesh, int eid, const std::vector<int>& node_ids);

    //! get the shape of the element
    Core::FE::CellType shape() override { return Core::FE::CellType::hex20; }

    //! compute local coordinates of the element for given global coordinates
    void local_coordinates(
        const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<3, 1>& rst) override;
  };

  /// hex27 element handle
  class Hex27ElementHandle : public QuadraticElementHandle
  {
   public:
    //! constructor
    Hex27ElementHandle(Mesh& mesh, int eid, const std::vector<int>& node_ids);

    //! get the shape of the element
    Core::FE::CellType shape() override { return Core::FE::CellType::hex27; }

    //! compute local coordinates of the element for given global coordinates
    void local_coordinates(
        const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<3, 1>& rst) override;
  };

  /// tet10 element handle
  class Tet10ElementHandle : public QuadraticElementHandle
  {
   public:
    //! constructor
    Tet10ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nids);

    //! get the shape of the element
    Core::FE::CellType shape() override { return Core::FE::CellType::tet10; }

    //! compute local coordinates of the element for given global coordinates
    void local_coordinates(
        const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<3, 1>& rst) override;
  };

  /// wedge15 element handle
  class Wedge15ElementHandle : public QuadraticElementHandle
  {
   public:
    //! constructor
    Wedge15ElementHandle(Mesh& mesh, int eid, const std::vector<int>& node_ids);

    //! get the shape of the element
    Core::FE::CellType shape() override { return Core::FE::CellType::wedge15; }

    //! compute local coordinates of the element for given global coordinates
    void local_coordinates(
        const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<3, 1>& rst) override;
  };

}  // namespace Cut


FOUR_C_NAMESPACE_CLOSE

#endif
