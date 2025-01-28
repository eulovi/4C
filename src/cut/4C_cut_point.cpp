// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_mesh.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_point_impl.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_tolerance.hpp"
#include "4C_global_data.hpp"

#include <string>

#ifdef CLN_CALC_OUTSIDE_KERNEL_POINT
#include "4C_cut_clnwrapper.hpp"
#endif

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Point* Cut::Point::new_point(
    Mesh& mesh, const double* x, double t, Edge* cut_edge, Side* cut_side, double tolerance)
{
  Point* p = mesh.new_point(x, cut_edge, cut_side, tolerance);
  p->position(Point::oncutsurface);
  // p->t( cut_edge, t );
  return p;
}

void Cut::Point::add_edge_intersection(Edge* first, Edge* second,
    const std::pair<Side*, Edge*>& original_cut_pair, const std::string& extra_msg)
{
#if CUT_CREATION_INFO
  std::stringstream msg;
  msg << extra_msg;
  msg << "//Added because of LATER edge-edge intersection of edge " << first->Id() << " and "
      << second->Id();
#endif

  // note: can be optimized
  const plain_side_set& first_sides = first->sides();
  for (plain_side_set::const_iterator j = first_sides.begin(); j != first_sides.end(); ++j)
  {
    Side* s = *j;
#if CUT_CREATION_INFO
    AddCreationInfo(std::make_pair(s, second), msg.str());
#endif
    add_pair(s, second, original_cut_pair);
  }
  const plain_side_set& second_sides = second->sides();
  for (plain_side_set::const_iterator j = second_sides.begin(); j != second_sides.end(); ++j)
  {
    Side* s = *j;
#if CUT_CREATION_INFO
    AddCreationInfo(std::make_pair(s, first), msg.str());
#endif
    add_pair(s, first, original_cut_pair);
  }
}

void Cut::Point::add_edge_intersection(Edge* first, Edge* second, Side* original_side,
    Edge* original_edge, const std::string& extra_msg)
{
  add_edge_intersection(first, second, std::make_pair(original_side, original_edge), extra_msg);
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Point* Cut::Point::insert_cut(Edge* cut_edge, Side* cut_side, Node* n)
{
  Point* p = n->point();
  const plain_edge_set& edges = n->edges();
  for (plain_edge_set::const_iterator i = edges.begin(); i != edges.end(); ++i)
  {
    Edge* e = *i;
    p->add_edge(e);
  }
  if (cut_edge != nullptr) p->add_edge(cut_edge);
  if (cut_side != nullptr) p->add_side(cut_side);

  if ((cut_side != nullptr) && (cut_edge != nullptr)) p->add_pair(cut_side, cut_edge);
  p->position(Point::oncutsurface);

#if CUT_CREATION_INFO
  std::pair<Side*, Edge*> cut_pair = std::make_pair(cut_side, cut_edge);
  p->AddMergedPair(cut_pair, nullptr);
#endif

  return p;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Point::Point(unsigned pid, double tolerance)
    : pid_(pid),
      position_(undecided),
      tol_(tolerance)
#if CUT_CREATION_INFO
      ,
      merged_to_(nullptr)
#endif
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::Point::add_edge(Edge* cut_edge)
{
  cut_edges_.insert(cut_edge);

  // reverse add
  cut_edge->add_point(this);

  const plain_side_set& edge_sides = cut_edge->sides();
  for (plain_side_set::const_iterator i = edge_sides.begin(); i != edge_sides.end(); ++i)
  {
    Side* s = *i;
    add_side(s);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::Point::add_side(Side* s)
{
  cut_sides_.insert(s);

  // revers add
  s->add_point(this);

  const plain_element_set& elements = s->elements();
  for (plain_element_set::const_iterator i = elements.begin(); i != elements.end(); ++i)
  {
    Element* e = *i;
    add_element(e);
  }
}

/*-----------------------------------------------------------------------------------*
 *    Identifies the edges that are cut by considered point and given point
 *-----------------------------------------------------------------------------------*/
void Cut::Point::common_edge(Point* other, plain_edge_set& edges)
{
  for (plain_edge_set::iterator i = cut_edges_.begin(); i != cut_edges_.end(); ++i)
  {
    Edge* e = *i;
    if (other->is_cut(e))
    {
      edges.insert(e);
    }
  }
}

/*-----------------------------------------------------------------------------------*
 *    Identifies the sides that are cut by considered point and given point
 *-----------------------------------------------------------------------------------*/
void Cut::Point::common_side(Point* other, plain_side_set& sides)
{
  for (plain_side_set::iterator i = cut_sides_.begin(); i != cut_sides_.end(); ++i)
  {
    Side* e = *i;
    if (other->is_cut(e))
    {
      sides.insert(e);
    }
  }
}

void Cut::Point::cut_edge(Side* side, Line* other_line, std::vector<Edge*>& matches)
{
  for (plain_edge_set::iterator i = cut_edges_.begin(); i != cut_edges_.end(); ++i)
  {
    Edge* e = *i;
    if (e->at_side(side))
    {
      if (e->begin_node()->point() == this or e->end_node()->point() == this)
      {
        if (not other_line->on_edge(e))
        {
          matches.push_back(e);
        }
      }
      else
      {
        matches.push_back(e);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Line* Cut::Point::cut_line(const Impl::PointLineFilter& filter, bool unique)
{
  Line* line_found = nullptr;
  for (plain_line_set::iterator i = lines_.begin(); i != lines_.end(); ++i)
  {
    Line* line = *i;
    if (filter(line))
    {
      if (line_found == nullptr)
      {
        line_found = line;
        if (not unique)
        {
          break;
        }
      }
      else
      {
        FOUR_C_THROW("not unique");
      }
    }
  }
  return line_found;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Line* Cut::Point::cut_line(Line* line, const Impl::PointLineFilter& filter, bool unique)
{
  Impl::ExcludeLineFilter f(line, filter);
  return cut_line(f, unique);
}

void Cut::Point::cut_lines(const Impl::PointLineFilter& filter, plain_line_set& cut_lines)
{
  for (plain_line_set::iterator i = lines_.begin(); i != lines_.end(); ++i)
  {
    Line* line = *i;
    if (filter(line))
    {
      cut_lines.insert(line);
    }
  }
}

void Cut::Point::cut_lines(Side* side, plain_line_set& cut_lines)
{
  Impl::SideCutFilter filter(side);
  Point::cut_lines(filter, cut_lines);
}


double Cut::Point::t(Edge* edge, const Core::LinAlg::Matrix<3, 1>& coord)
{
  std::map<Edge*, double>::iterator i = t_.find(edge);
  if (i == t_.end())
  {
    Core::LinAlg::Matrix<3, 1> x(coord.data());

    Point* p1 = edge->begin_node()->point();
    Point* p2 = edge->end_node()->point();

    if (p1 == p2) return 0;

    if (p1 == this) return -1;

    if (p2 == this) return 1;

    Core::LinAlg::Matrix<3, 1> x1;
    Core::LinAlg::Matrix<3, 1> x2;

    p1->coordinates(x1.data());
    p2->coordinates(x2.data());


#ifdef CLN_CALC_OUTSIDE_KERNEL_POINT

    Core::LinAlg::Matrix<3, 1, Core::CLN::ClnWrapper> _cln;
    int prec_old = Core::CLN::ClnWrapper::precision_;
    Core::CLN::ClnWrapper::precision_ = 50;
    int prec = Core::CLN::ClnWrapper::GetPrecision();
    ConvDoubleCLN(x, x_cln., prec);

    Core::LinAlg::Matrix<3, 1, Core::CLN::ClnWrapper> x1_cln;
    Core::LinAlg::Matrix<3, 1, Core::CLN::ClnWrapper> x2_cln;

    ConvDoubleCLN(x1, x1_cln., prec);
    ConvDoubleCLN(x2, x2_cln., prec);

    x_cln.update(-1.0, x1_cln, 1.0);
    x2_cln.update(-1.0, x1_cln, 1.0);

    Core::CLN::ClnWrapper l1_cln = x_cln.norm2();
    Core::CLN::ClnWrapper l2_cln = x2_cln.norm2();

    if (Core::MathOperations<Core::CLN::ClnWrapper>::abs(l2_cln) <
        (p1->Tolerance() + p2->Tolerance()))
    {
      FOUR_C_THROW("edge with no length");
    }


    Core::CLN::ClnWrapper z_cln = l1_cln / l2_cln;

    x_cln.update(-z_cln, x2_cln, 1.0);

    if (x_cln.norm2() > (Tolerance() + p1->Tolerance() + p2->Tolerance()))
    {
      std::stringstream str;
      str << "point id " << this->Id() << "\n";
      str << "point not on edge, no edge position: " << x.norm2()
          << " (Tol = " << Tolerance() + p1->Tolerance() + p2->Tolerance() << ")"
          << "\n"
          << x << x1 << x2 << "\n";
      FOUR_C_THROW(str.str());
    }
    // transformation to the parameter space coordinate t of the edge (between -1 and 1)
    Core::CLN::ClnWrapper t_cln = 2.0 * z_cln - 1.0;
    double t = Core::MathOperations<Core::CLN::ClnWrapper>::GetDouble(t_cln);
    Core::CLN::ClnWrapper::precision_ = prec_old;  // restoring previous cln precision

#else
    x.update(-1, x1, 1);
    x2.update(-1, x1, 1);

    double l1 = x.norm2();
    double l2 = x2.norm2();

    if (fabs(l2) < (p1->tolerance() + p2->tolerance()))
    {
      FOUR_C_THROW("edge with no length");
    }

    double z = l1 / l2;

    x.update(-z, x2, 1);
    // one could think of choosing a tighter tolerance here, but why?
    if (x.norm2() > (tolerance() + p1->tolerance() + p2->tolerance()))
    {
      std::ofstream file("point_not_on_on_edge.pos");
      Cut::Output::gmsh_edge_dump(file, edge, std::string("Edge"));
      Cut::Output::gmsh_point_dump(file, this, this->id(), std::string("Point"), false, nullptr);
      file.close();

      std::stringstream str;
      str << "point id " << this->id() << "\n";
      str << "point not on edge, no edge position: " << x.norm2()
          << " (Tol = " << tolerance() + p1->tolerance() + p2->tolerance() << ")"
          << "\n"
          << x << x1 << x2 << "\n";
      FOUR_C_THROW(str.str());
    }
    // transformation to the parameter space coordinate t of the edge (between -1 and 1)
    double t = 2.0 * z - 1.0;
#endif
    const PointPositionSet& current_cut_points = edge->cut_points();
    // try to find if point with the same local coordinate already exist on the edge and correct its
    // location
    for (PointPositionSet::const_iterator it = current_cut_points.begin();
        it != current_cut_points.end(); ++it)
    {
      // just query, since they are already added
      if ((*it)->t(edge) == t)
      {
        // if the local coordinate of the point corresponds to the end node of the edge (possible
        // due to precision), we shift endpoint a bit outside to preserve correct order and make
        // them distinguishable
        if ((*it)->nodal_point(edge->nodes()) and ((t == 1.0 || t == -1.0)))
        {
          if (edge->begin_node()->point() == (*it))
            (*it)->t_[edge] = t - END_NODE_SHIFT_DISTANCE;
          else
            (*it)->t_[edge] = t + END_NODE_SHIFT_DISTANCE;
#if EXTENDED_CUT_DEBUG_OUTPUT
          std::cout << "Inserting point to the edge, close to its end points but without merging"
                    << std::endl;
#endif
        }
        else
        {
          // crash
          std::stringstream str;
          str << "Two points with the same local coordinates on the edge. " << id()
              << "will not be added to a list because of existing " << (*it)->id()
              << " in pointposition list\n";
          str << "Distance between points is " << Cut::distance_between_points(this, *it)
              << std::endl;
          str << this << std::endl;
          str << *it << std::endl;
          str << "The local coordinates are " << std::setprecision(15) << (*it)->t(edge) << "and"
              << t << std::endl;

          if ((*it)->nodal_point(edge->nodes()))
            str << "Shifting edge point second time\n";
          else
            str << "Two non end-point have same coordinates\n";

          std::ofstream file("t_failed.pos");
          Cut::Output::gmsh_edge_dump(file, edge);
          file.close();
          dump_connectivity_info();
          (*it)->dump_connectivity_info();
          FOUR_C_THROW(str.str());
        }
      }
    }
    t_[edge] = t;
    return t;
  }
  return i->second;
}

double Cut::Point::t(Edge* edge)
{
  Core::LinAlg::Matrix<3, 1> x;
  coordinates(x.data());
  double rv = Cut::Point::t(edge, x);
  return rv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::Point::intersection(plain_edge_set& edges)
{
  plain_edge_set intersection;
  std::set_intersection(cut_edges_.begin(), cut_edges_.end(), edges.begin(), edges.end(),
      std::inserter(intersection, intersection.begin()));
  std::swap(edges, intersection);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::Point::intersection(plain_side_set& sides)
{
  plain_side_set intersection;
  std::set_intersection(cut_sides_.begin(), cut_sides_.end(), sides.begin(), sides.end(),
      std::inserter(intersection, intersection.begin()));
  std::swap(sides, intersection);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::Point::intersection(plain_facet_set& facets)
{
  plain_facet_set intersection;
  std::set_intersection(facets_.begin(), facets_.end(), facets.begin(), facets.end(),
      std::inserter(intersection, intersection.begin()));
  std::swap(facets, intersection);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::Point::intersection(plain_element_set& elements)
{
  plain_element_set intersection;
  std::set_intersection(cut_elements_.begin(), cut_elements_.end(), elements.begin(),
      elements.end(), std::inserter(intersection, intersection.begin()));
  std::swap(elements, intersection);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::Point::nodal_point(const std::vector<Node*>& nodes) const
{
  for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    Node* n = *i;
    if (n->point() == this)
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::Point::almost_nodal_point(const std::vector<Node*>& nodes, double tolerance) const
{
  Core::LinAlg::Matrix<3, 1> np;
  Core::LinAlg::Matrix<3, 1> p;

  this->coordinates(p.data());

  for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    Node* n = *i;
    n->point()->coordinates(np.data());
    np.update(-1, p, 1);

    if (np.norm2() <= tolerance)
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Node* Cut::Point::cut_node()
{
  for (plain_edge_set::iterator i = cut_edges_.begin(); i != cut_edges_.end(); ++i)
  {
    Edge* e = *i;
    const std::vector<Node*>& nodes = e->nodes();
    for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
    {
      Node* n = *i;
      if (n->point() == this)
      {
        return n;
      }
    }
  }
  return nullptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::Point::position(Point::PointPosition pos)
{
  if (position_ != pos)
  {
    //  safety check, if the position of a facet changes from one side to the other
    FOUR_C_ASSERT(is_cut_position_unchanged(position_, pos),
        "Are you sure that you want to change the point-position from inside to outside or vice "
        "versa?");

    // do not overwrite oncutsurface points
    if (position_ == Point::oncutsurface) return;

    // change position for points just in case of undecided point and do not change oncutsurface
    // points
    if (position_ == undecided)
    {
      position_ = pos;
      if (pos == Point::outside or pos == Point::inside)
      {
        for (plain_facet_set::iterator i = facets_.begin(); i != facets_.end(); ++i)
        {
          Facet* f = *i;
          if (f != nullptr) f->position(pos);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::Point::has_associated_boundary_cell_facet()
{
  for (plain_facet_set::const_iterator ifacet = facets_.begin(); ifacet != facets_.end(); ++ifacet)
  {
    if ((*ifacet)->on_boundary_cell_side()) return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Side* Cut::Point::cut_side(Side* side, Point* other)
{
  Side* found_side = nullptr;
  for (plain_side_set::iterator i = cut_sides_.begin(); i != cut_sides_.end(); ++i)
  {
    Side* s = *i;
    if (s != side and other->is_cut(s))
    {
      if (found_side == nullptr)
      {
        found_side = s;
      }
      else
      {
        FOUR_C_THROW("side not unique");
      }
    }
  }
  return found_side;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Cut::Point> Cut::create_point(
    unsigned pid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance)
{
  const PointFactory factory;
  const int probdim = Global::Problem::instance()->n_dim();
  return factory.create_point(pid, x, cut_edge, cut_side, tolerance, probdim);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::Impl::SideCutFilter::operator()(Line* line) const { return line->is_internal_cut(side_); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::Impl::SideElementCutFilter::operator()(Line* line) const
{
  // return line->IsCut( side_ ) and line->IsCut( element_ );
  return line->is_internal_cut(side_) and line->is_cut(element_);
}

bool Cut::Point::is_cut(Side* side, Edge* edge)
{
  bool ret = cut_pairs_.find(std::pair<Side*, Edge*>(side, edge)) != cut_pairs_.end();
  return ret;
}

void Cut::Point::add_pair(Side* side, Edge* edge, const std::pair<Side*, Edge*>& original_pair)
{
  cut_pairs_.insert(std::pair<Side*, Edge*>(side, edge));
  // does not overwrite existing value
#if CUT_CREATION_INFO
  cut_pairs_info_.insert(
      std::make_pair((std::pair<Side*, Edge*>(side, edge)), std::make_pair(original_pair, this)));
#endif
}

void Cut::Point::add_pair(const std::pair<Side*, Edge*>& pair,
    const std::pair<Side*, Edge*>& original_pair, Point* source)
{
  cut_pairs_.insert(pair);
#if CUT_CREATION_INFO
  cut_pairs_info_.insert(std::make_pair(pair, std::make_pair(original_pair, source)));
#endif
}


void Cut::Point::add_pair(Side* side, Edge* edge)
{
  add_pair(side, edge, std::make_pair(side, edge));
}

bool Cut::Point::is_cut(Side* s1, Side* s2)
{
  const std::vector<Edge*>& edges_s1 = s1->edges();
  const std::vector<Edge*>& edges_s2 = s2->edges();

  for (std::vector<Edge*>::const_iterator it = edges_s2.begin(); it != edges_s2.end(); ++it)
  {
    if (this->is_cut(s1, (*it))) return true;
  }

  for (std::vector<Edge*>::const_iterator it = edges_s1.begin(); it != edges_s1.end(); ++it)
  {
    if (this->is_cut(s2, (*it))) return true;
  }

  // no combination intersects
  return false;
}

void Cut::Point::dump_connectivity_info()
{
  std::stringstream p;
  p << "Point" << id() << "ConnectivityDump.pos";
  std::string fname = p.str();
  std::ofstream file(fname.c_str());

#if (not CUT_CREATION_INFO)

  file << "// Enable CREATION_INFO flag in the src/cut/cut_tolerance.H to get full output of"
          "how this point was created"
       << std::endl;
#else

  int counter = 0;
  plain_edge_set part_edges;
  plain_side_set part_sides;
  std::pair<Side*, Edge*> or_cut;  // original intersection

  file << "//==========================================\n"
          "//== Additional cut pair independent info ==\n"
          "//=========================================="
       << std::endl;

  file << GetAdditionalCreationInfo() << std::endl;

  file << "// This point coordinates" << std::endl;
  std::stringstream point_name;
  point_name << "Point" << Id();
  Cut::Output::GmshPointDump(file, this, Id(), point_name.str(), false, nullptr);

  if (merged_to_)
  {
    file << "This point was merged into" << merged_to_->Id() << "\n";
    file.close();
  }
  else
  {
    for (std::set<std::pair<Side*, Edge*>>::iterator it = cut_pairs_.begin();
        it != cut_pairs_.end(); ++it, ++counter)
    {
      std::stringstream prefix;
      Point* p_from = cut_pairs_info_[*it].second;
      if (p_from != this)
      {
        prefix << "_Merged_";
        Point* p_next = p_from;
        do
        {
          prefix << p_next->Id();
          if (p_next != this) prefix << "->";
          p_from = p_next;
          p_next = p_next->merged_to_;
        } while (p_next != nullptr);


        prefix << "_";

        std::vector<Point*>::iterator pt =
            std::find(real_merged_points_.begin(), real_merged_points_.end(), p_from);
        if (pt == real_merged_points_.end()) FOUR_C_THROW("Not possible");
      }

      std::string pre = prefix.str();
      file << "\n// ===Start new pair info===";
      if (IsReallyCut(*it))
      {
        file << "// was really cut by this combination\n";
        or_cut = *it;
      }
      else
      {
        const std::pair<Side*, Edge*>& original_connection = this->AddedFrom(*it);
        file << "// original intersection was and because of it" << std::endl;
        std::stringstream side_section_name;
        side_section_name << "OriginalSide" << pre << counter;
        Cut::Output::gmsh_new_section(file, side_section_name.str());
        Cut::Output::gmsh_side_dump(file, original_connection.first, false, nullptr);
        Cut::Output::gmsh_end_section(file, false);

        std::stringstream edge_section_name;
        edge_section_name << "OriginalEdge" << pre << counter;
        Cut::Output::gmsh_new_section(file, edge_section_name.str());
        Cut::Output::GmshEdgeDump(file, original_connection.second, false, nullptr);
        Cut::Output::gmsh_end_section(file, false);
        or_cut = original_connection;
        file << "// " << GetCreationInfo(or_cut) << "\n";
        file << "// following connection is added later" << std::endl;
      }
      std::stringstream edge_section_name;
      edge_section_name << "Edge" << pre << counter;
      Cut::Output::gmsh_new_section(file, edge_section_name.str());
      Cut::Output::GmshEdgeDump(file, (*it).second);
      Cut::Output::gmsh_end_section(file, false);
      std::stringstream side_section_name;
      side_section_name << "Side" << pre << counter;
      Cut::Output::gmsh_new_section(file, side_section_name.str());
      Cut::Output::gmsh_side_dump(file, (*it).first);
      Cut::Output::gmsh_end_section(file, false);
      file << "// " << GetCreationInfo(*it) << "\n";

      std::map<std::pair<Side*, Edge*>, std::pair<Core::LinAlg::Matrix<3, 1>, bool>>::iterator
          merged = merged_points_info_.find(or_cut);
      if (merged != merged_points_info_.end())
      {
        std::pair<Core::LinAlg::Matrix<3, 1>, bool>& merged_coord = merged->second;
        if (merged_coord.second)
        {
          file << "// original point coordinates from real intersection \n";
          file << "View \"MergedPoint" << Id() << "_" << counter << "\"{\n";
          file << "SP (";
          Cut::Output::GmshWriteCoords(file, merged_coord.first, false, nullptr);
          file << "){";
          file << 0;
          Cut::Output::gmsh_end_section(file, false);
          Cut::Output::gmsh_end_section(file, false);
          file << "// Norm2 between point and merged is "
               << Cut::DistanceBetweenPoints(this, merged_coord.first) << std::endl;
        }
        else
          file << "// Another point was merged here during insert cut\n";
      }
      else
      {
        file << "// This was result of pointpool merge, but this was due some topological "
                "connection added before!";
      }
      part_sides.insert((*it).first);
      part_edges.insert((*it).second);
      file << "\n";
    }

    for (plain_edge_set::iterator it = cut_edges_.begin(); it != cut_edges_.end(); ++it, ++counter)
    {
      if (part_edges.find(*it) == part_edges.end())
      {
        std::stringstream section_name;
        section_name << "IndividualEdge" << counter;
        Cut::Output::gmsh_new_section(file, section_name.str());
        Cut::Output::GmshEdgeDump(file, (*it));
        Cut::Output::gmsh_end_section(file, false);
        if (not NodalPoint((*it)->Nodes())) file << "// it is not a nodal point!" << std::endl;
      }
    }

    for (plain_side_set::iterator it = cut_sides_.begin(); it != cut_sides_.end(); ++it, ++counter)
    {
      if (part_sides.find(*it) == part_sides.end())
      {
        std::stringstream section_name;
        section_name << "IndividualSide" << counter;
        Cut::Output::gmsh_new_section(file, section_name.str());
        Cut::Output::gmsh_side_dump(file, (*it));
        Cut::Output::gmsh_end_section(file, false);
        if (not NodalPoint((*it)->Nodes())) file << " // it is not a nodal point!" << std::endl;
      }
    }
  }
#endif
  file.close();
}

/// Remove info about this point from all the cut_sides and cut_edges
void Cut::Point::remove_connectivity_info()
{
  for (plain_side_set::iterator it = cut_sides_.begin(); it != cut_sides_.end(); ++it)
    (*it)->remove_point(this);
  for (plain_edge_set::iterator it = cut_edges_.begin(); it != cut_edges_.end(); ++it)
    (*it)->remove_point(this);
  /* Dont remove from cut_element, because it does not store cut_points inside itself */

  cut_pairs_.clear();
  cut_sides_.clear();
  cut_edges_.clear();
  cut_elements_.clear();
}

void Cut::Point::replace(Point* p)
{
  p->cut_pairs_.insert(cut_pairs_.begin(), cut_pairs_.end());
  p->cut_sides_.insert(cut_sides_.begin(), cut_sides_.end());
  p->cut_edges_.insert(cut_edges_.begin(), cut_edges_.end());
  p->cut_elements_.insert(cut_elements_.begin(), cut_elements_.end());

  for (Side* s : cut_sides_)
  {
    s->remove_point(this);
    s->add_point(p);
    std::set<std::set<Point*>>& parallel_cut_surfaces = s->get_parallel_cut_surfaces();
    std::set<std::set<Point*>> new_surfaces;
    std::transform(parallel_cut_surfaces.begin(), parallel_cut_surfaces.end(),
        std::inserter(new_surfaces, new_surfaces.begin()),
        [&](std::set<Point*> surface)
        {
          surface.erase(this);
          surface.insert(p);
          return surface;
        });
    parallel_cut_surfaces.swap(new_surfaces);
  }

  Core::LinAlg::Matrix<3, 1> real_coord;
  coordinates(real_coord.data());
  for (Edge* e : cut_edges_)
  {
    e->remove_point(this);
    // preserve location on this edge with respect to other points
    p->t(e, real_coord);
    e->add_point(p);
  }

  cut_pairs_.clear();
  cut_sides_.clear();
  cut_edges_.clear();
  cut_elements_.clear();

  plain_line_set replaces;
  for (plain_line_set::iterator it = lines_.begin(); it != lines_.end();)
  {
    Line* existing_line = *it;
    Point* other = existing_line->other_point(this);
    Line* replace_line = other->common_line(p);
    // no such line exist, modify current to include new endpoint
    if (replace_line == nullptr)
    {
      existing_line->replace(this, p);
      ++it;
    }
    // use already existed line as replacement
    else
    {
      const plain_side_set& lines_sides = existing_line->cut_sides();
      for (auto& s : lines_sides)
      {
        auto& side_lines = const_cast<std::vector<Line*>&>(s->cut_lines());
        auto erase_it = std::find(side_lines.begin(), side_lines.end(), existing_line);
        if (erase_it != side_lines.end())
        {
          side_lines.erase(erase_it);
          if (std::find(side_lines.begin(), side_lines.end(), replace_line) == side_lines.end())
            side_lines.push_back(replace_line);
        }
      }
      replaces.insert(replace_line);
#ifdef CUT_USE_SORTED_VECTOR
      it = lines_.ierase(it);
#else
      it = lines_.erase(it);
#endif
    }
  }

  p->lines_.insert(lines_.begin(), lines_.end());
  p->lines_.insert(replaces.begin(), replaces.end());
  lines_.clear();
  p->merged_points_.push_back(this);
}

void Cut::Point::merge(Point* dest)
{
#if CUT_CREATION_INFO
  merged_to_ = dest;
  dest->real_merged_points_.push_back(dest);
#endif
}

void Cut::Point::remove_side(Side* side)
{
  side->remove_point(this);
  cut_sides_.erase(side);
}

void Cut::Point::remove_edge(Edge* edge)
{
  edge->remove_point(this);
  cut_edges_.erase(edge);
}

void Cut::Point::erased_containing_cut_pairs(Side* side)
{
  for (std::set<std::pair<Side*, Edge*>>::iterator i = cut_pairs_.begin(); i != cut_pairs_.end();)
  {
    std::pair<Side*, Edge*> cut_pair = *i;
    if (cut_pair.first == side)
    {
      i = cut_pairs_.erase(i);
    }
    else
    {
      ++i;
    }
  }
}

void Cut::Point::erased_containing_cut_pairs(Edge* edge)
{
  for (std::set<std::pair<Side*, Edge*>>::iterator i = cut_pairs_.begin(); i != cut_pairs_.end();)
  {
    std::pair<Side*, Edge*> cut_pair = *i;
    if (cut_pair.second == edge)
    {
      i = cut_pairs_.erase(i);
    }
    else
    {
      ++i;
    }
  }
}


#if CUT_CREATION_INFO

void Cut::Point::AddCreationInfo(const std::pair<Side*, Edge*>& cut_pair, const std::string& info)
{
  typedef std::map<std::pair<Side*, Edge*>, std::string>::iterator info_iterator;
  std::pair<info_iterator, bool> inserted = creation_info_.insert(std::make_pair(cut_pair, info));
  if (not inserted.second)
  {
    inserted.first->second.append("// AND ").append(info);
  }
}

void Cut::Point::AddAdditionalCreationInfo(const std::string& info)
{
  additional_creation_info_.append("// And more: \n");
  additional_creation_info_.append(info);
}

const std::string& Cut::Point::GetAdditionalCreationInfo() { return additional_creation_info_; }


void Cut::Point::AddMergedPair(
    const std::pair<Side*, Edge*>& inter, Core::LinAlg::Matrix<3, 1>* coord)
{
  std::map<std::pair<Side*, Edge*>, std::pair<Core::LinAlg::Matrix<3, 1>, bool>>::iterator merged =
      merged_points_info_.find(inter);
  if (merged != merged_points_info_.end())
  {
    // actually this can happen
    // FOUR_C_THROW("This should not happen"); //
  }
  else if (coord)
    merged_points_info_[inter] = std::make_pair(*coord, true);
  else
  {
    Core::LinAlg::Matrix<3, 1> dummy;
    merged_points_info_[inter] = std::make_pair(dummy, false);
  }
}

#endif

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::Impl::SideSideCutFilter::operator()(Line* line) const
{
  return line->is_cut(side1_) and line->is_cut(side2_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::Edge* Cut::Point::common_cut_edge(Side* side)
{
  for (plain_edge_set::iterator i = cut_edges_.begin(); i != cut_edges_.end(); ++i)
  {
    Edge* e = *i;
    if (e->is_cut(side))
    {
      return e;
    }
  }
  return nullptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& stream, Cut::Point& point)
{
  point.print(stream);
  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& stream, Cut::Point* point)
{
  point->print(stream);
  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::find_common_elements(const std::vector<Point*>& element, plain_element_set& elements)
{
  std::vector<Point*>::const_iterator ie = element.begin();
  elements = (*ie)->elements();
  for (++ie; ie != element.end(); ++ie)
  {
    Point* p = *ie;
    p->intersection(elements);
    if (elements.size() == 0)
    {
      break;
    }
  }
}

double Cut::distance_between_points(Point* p1, Point* p2)
{
  Core::LinAlg::Matrix<3, 1> p1_x;
  p1->coordinates(p1_x.data());
  Core::LinAlg::Matrix<3, 1> p2_x;
  p2->coordinates(p2_x.data());
  return distance_between_points(p1_x, p2_x);
}

double Cut::distance_between_points(Point* p1, const Core::LinAlg::Matrix<3, 1>& coord_b)
{
  Core::LinAlg::Matrix<3, 1> p1_x;
  p1->coordinates(p1_x.data());
  return distance_between_points(p1_x, coord_b);
}


bool Cut::is_cut_position_unchanged(Point::PointPosition position, Point::PointPosition pos)
{
  if ((position == Point::inside and pos == Point::outside) or
      (position == Point::outside and pos == Point::inside))
    return false;
  else
    return true;
}

template <unsigned prob_dim>
void Cut::ConcretePoint<prob_dim>::move_point(const double* new_coord)
{
  Core::LinAlg::Matrix<3, 1> coordm(new_coord, true);
  if (Cut::distance_between_points(this, coordm) > BOXOVERLAP)
    FOUR_C_THROW("The point is not allowed to be moved that much");

#if CUT_CREATION_INFO
  std::stringstream info;
  info << "// This point was moved to the new coordinates " << coordm << std::endl;
  AddAdditionalCreationInfo(info.str());
#endif
  std::copy(new_coord, new_coord + prob_dim, this->x_);
}

FOUR_C_NAMESPACE_CLOSE
