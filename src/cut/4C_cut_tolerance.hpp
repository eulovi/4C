// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_TOLERANCE_HPP
#define FOUR_C_CUT_TOLERANCE_HPP

#include "4C_config.hpp"

#define BASICTOL 1e-14  // basic tolerance for the cut, limited by double precicion!

// derived tolerances
//--> do not change the relation of these tolerances from the basictolerance, if there isn't a
// reasonable argument to do that!!!

// to check l2-norm of vector with 3-entries of magnitude 1.0!
#define LINSOLVETOL (BASICTOL / sqrt(3) * 10.0)
// sqrt(3) because of l2 norm with 3 values in vector
// 10 because scaled with norminf(x) || if norminf(x)=1.0 --> for double precicion 1.(15 digits)
// REMARK: using this tolerance for increments in reference coordinate system norminf(x)=1!!!

// to define general tolerance for position represented by a vector of doubles!
#define POSITIONTOL (BASICTOL * 10.0)


#define REFERENCETOL \
  (BASICTOL * 100.0)  // tolerance to check coordinates in the relative coordinate system (assumtion
                     // is that abs(max(value)) = 1.0)!

// Not in use anymore.
#define VOLUMETOL \
  1e-14  // JUST FOR TESSELATION: the tolerance when created integration cells are discarded, the
         // boundary cell however will be created

// special tolerances for some special cut situations --> we should define them relative to our
// BASICTOL!!!

// tolerance to identify points of volumecells between processors
#define PARALLEL_COORD_TOL 1e-10

// Special tolerance to distinguish planar facets from curved ones.
#define PLANARTOL 1e-6

#define BOXOVERLAP 1e-7

// if cross product constructed from any 3 points are less than this,
// then these points are on the same line
#define TOL_POINTS_ON_LINE 1e-14

#define TOL_POINTS_ON_LINE_FOR_DELETING 1e-11

// In equation of plane ax+by+cz=d for a facet, this tolerance is used to check
// whether any of the four coefficients are zero
#define TOL_EQN_PLANE 1e-12

// For moment fitting, during line integration of base functions, if normal componenet of
// lines of less than this, they are considered zero
#define TOL_LINE_NORMAL 1e-12

// Obtaining  internal Gauss points in directdivergence involves division by
// coefficent 'a' in equation of plane. A facet whose 'a' is more than this tolerance is
// chosen as reference plane to avoid a large number (and hence large error)
#define REF_PLANE_DIRDIV 1e-3

// We may get volumecells with all facets on a plane
// in such case directdivergence may produce negative volume, because it works only for closed
// geometries if volume of cell is less than this, we do not generate gauss rules
#define REF_VOL_DIRDIV 1e-6

// Tolerance use in the Selfcut for detection of the position insider or outside
#define SELF_CUT_POS_TOL 1e-8

// If the area of boundarycell is less than this value, then it is eliminated from boundary
// integration Such small area can lead to zero or a very small negative value in the boundary cell
// mapping's Jacobian This is only for the case of Direct Divergence!
#define REF_AREA_BCELL 1e-13  // Make it viable for GLOBAL!!!

// If the distance of the line between the main Gauss points and its projected point is less than
// this value, then this main Gauss point is not considered
#define JAC_LINE_TOL 1e-14

// Defines the order of the Gaussrule used for the direct divergence.
//  To make it fit also the DivergenceCells have to be fully triangulated for the mesh as is done in
//  the LS. See the DirectDiv algorithm in facetintegration.
#define DIRECTDIV_GAUSSRULE 7  // 7 for test cases

// Construct the integration rules directly on the local coordinate system of the element
// #define LOCAL

// utilized quadrature compression (for Details see paper by Sudhakar et al.)
// #define QUADCOMP

// if we want to use cln outside of cut kernel (e.g. for calculation of local position of point on
// edge or arc length in the pointgraph). Should be ajusted for memory containter in order to be
// usable
// #define CLN_CALC_OUTSIDE_KERNEL


// limiting error, until which precision will be increased
#define CLN_LIMIT_ERROR 1e-20

// limiting error after which double will be switched to cln
#define DOUBLE_LIMIT_ERROR (4 * 1e-15)

// whether we run on double + (soemtimes) cln or double + (always) cln
#define DOUBLE_PLUS_CLN_COMPUTE true

// global tolerance for detecting sides near the point in the cut_kernel
#define SIDE_DETECTION_TOLERANCE 1e-14

// global tolerance for detecting sides near the point in the cut_kernel
#define INSIDE_OUTSIDE_TOLERANCE 0.0

// new general mergin tolerance
#define MERGING_TOLERANCE 1e-14

// tolerance for merging non-topologically connected points
#define NON_TOPOLOGICAL_TOLERANCE 1e-14

// tolerance for mergin topologically connected points
#define TOPOLOGICAL_TOLERANCE 1e-14

// tolerance for detection if the edges are close enough to each other
#define EDGE_EDGE_TOLERANCE 1e-14

// minimum vaue to set up, when error of double calculation is equal 0.0
#define MINIMUM_DOUBLE_ERROR 1e-16

// scale of the tolerance, during initial loading of the mesh
#define NODAL_POINT_TOLERANCE_SCALE 5.0

// scale of the tolerance, during creation of points in the selcut
#define NODAL_POINT_TOLERANCE_SELFCUT_SCALE 2.5

// for direct_divergence.cpp volume of bounding box of facets with x-normal zero
#define MINIMUM_VOLUME_BB_FACETS 1e-13

// maximum area of the boundary cell, that might be negative due to numerical precision
#define AREA_BC_TOL 1e-12

// allowed error if volume tests
#define VOLUME_ERROR_TOL 1e-2

// distance with with to shift local position of the end points on the edge if it is equal to the
// othe cut_point local position
#define END_NODE_SHIFT_DISTANCE 1e-15

// cln memory manager

// #define CUSTOM_MEMORY_ALLOCATOR

#define DEBUG_MEMORY_ALLOCATION false

// define flags for effective cut development

// version for development, settign this to true will set all the needed constants also
#define CUT_DEVELOP false

#if CUT_DEVELOP
// all the basicl, useful cut output output, that is needed for development
#define EXTENDED_CUT_DEBUG_OUTPUT true
// whether to record information about how this point was created, its topology infromation,etc
#define CUT_CREATION_INFO true
#else
#define EXTENDED_CUT_DEBUG_OUTPUT false
#define CUT_CREATION_INFO false

#endif

#endif
