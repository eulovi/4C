// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_meshreader.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element_definition.hpp"
#include "4C_fem_general_fiber_node.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_io_domainreader.hpp"
#include "4C_io_elementreader.hpp"
#include "4C_io_exodus.hpp"
#include "4C_io_input_file.hpp"
#include "4C_io_value_parser.hpp"
#include "4C_rebalance.hpp"
#include "4C_rebalance_graph_based.hpp"
#include "4C_rebalance_print.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <utility>

FOUR_C_NAMESPACE_OPEN


namespace
{
  std::vector<std::shared_ptr<Core::FE::Discretization>> find_dis_node(
      const std::vector<Core::IO::ElementReader>& element_readers, int global_node_id)
  {
    std::vector<std::shared_ptr<Core::FE::Discretization>> list_of_discretizations;
    for (const auto& element_reader : element_readers)
      if (element_reader.has_node(global_node_id))
        list_of_discretizations.emplace_back(element_reader.get_dis());

    return list_of_discretizations;
  }

  void read_nodes(Core::IO::InputFile& input, const std::string& node_section_name,
      std::vector<Core::IO::ElementReader>& element_readers, int& max_node_id)
  {
    const int myrank = Core::Communication::my_mpi_rank(input.get_comm());
    if (myrank > 0) return;

    int line_count = 0;
    for (const auto& node_line : input.in_section_rank_0_only(node_section_name))
    {
      Core::IO::ValueParser parser{
          node_line.get_as_dat_style_string(), {.user_scope_message = "While reading node data: "}};
      auto type = parser.read<std::string>();

      if (type == "NODE")
      {
        int nodeid = parser.read<int>() - 1;
        parser.consume("COORD");
        auto coords = parser.read<std::vector<double>>(3);

        max_node_id = std::max(max_node_id, nodeid) + 1;
        std::vector<std::shared_ptr<Core::FE::Discretization>> dis =
            find_dis_node(element_readers, nodeid);

        for (const auto& di : dis)
        {
          // create node and add to discretization
          std::shared_ptr<Core::Nodes::Node> node =
              std::make_shared<Core::Nodes::Node>(nodeid, coords, myrank);
          di->add_node(node);
        }
      }
      // this node is a Nurbs control point
      else if (type == "CP")
      {
        int cpid = parser.read<int>() - 1;
        parser.consume("COORD");
        auto coords = parser.read<std::vector<double>>(3);
        double weight = parser.read<double>();

        max_node_id = std::max(max_node_id, cpid) + 1;
        if (cpid != line_count)
          FOUR_C_THROW(
              "Reading of control points {} failed: They must be numbered consecutive!!", cpid);
        std::vector<std::shared_ptr<Core::FE::Discretization>> diss =
            find_dis_node(element_readers, cpid);

        for (auto& dis : diss)
        {
          // create node/control point and add to discretization
          std::shared_ptr<Core::FE::Nurbs::ControlPoint> node =
              std::make_shared<Core::FE::Nurbs::ControlPoint>(cpid, coords, weight, myrank);
          dis->add_node(node);
        }
      }
      // this is a special node with additional fiber information
      else if (type == "FNODE")
      {
        enum class FiberType
        {
          Unknown,
          Angle,
          Fiber,
          CosyDirection
        };

        // read fiber node
        std::map<Core::Nodes::CoordinateSystemDirection, std::array<double, 3>> cosyDirections;
        std::vector<std::array<double, 3>> fibers;
        std::map<Core::Nodes::AngleType, double> angles;

        int nodeid = parser.read<int>() - 1;
        parser.consume("COORD");
        auto coords = parser.read<std::vector<double>>(3);
        max_node_id = std::max(max_node_id, nodeid) + 1;

        while (!parser.at_end())
        {
          auto next = parser.read<std::string>();

          if (next == "FIBER" + std::to_string(1 + fibers.size()))
          {
            fibers.emplace_back(parser.read<std::array<double, 3>>());
          }
          else if (next == "CIR")
          {
            cosyDirections[Core::Nodes::CoordinateSystemDirection::Circular] =
                parser.read<std::array<double, 3>>();
          }
          else if (next == "TAN")
          {
            cosyDirections[Core::Nodes::CoordinateSystemDirection::Tangential] =
                parser.read<std::array<double, 3>>();
          }
          else if (next == "RAD")
          {
            cosyDirections[Core::Nodes::CoordinateSystemDirection::Radial] =
                parser.read<std::array<double, 3>>();
          }
          else if (next == "HELIX")
          {
            angles[Core::Nodes::AngleType::Helix] = parser.read<double>();
          }
          else if (next == "TRANS")
          {
            angles[Core::Nodes::AngleType::Transverse] = parser.read<double>();
          }
        }

        // add fiber information to node
        std::vector<std::shared_ptr<Core::FE::Discretization>> discretizations =
            find_dis_node(element_readers, nodeid);
        for (auto& dis : discretizations)
        {
          auto node = std::make_shared<Core::Nodes::FiberNode>(
              nodeid, coords, cosyDirections, fibers, angles, myrank);
          dis->add_node(node);
        }
      }
      else
        FOUR_C_THROW("Unknown node type '{}'", type);

      ++line_count;
    }
  }
}  // namespace

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::MeshReader::MeshReader(
    Core::IO::InputFile& input, std::string node_section_name, MeshReaderParameters parameters)
    : comm_(input.get_comm()),
      input_(input),
      node_section_name_(std::move(node_section_name)),
      parameters_(std::move(parameters))
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::add_advanced_reader(
    std::shared_ptr<Core::FE::Discretization> dis, const std::string& sectionname)
{
  target_discretizations_.emplace_back(sectionname, dis);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::read_and_partition()
{
  // We need to track the max global node ID to offset node numbering and for sanity checks
  int max_node_id = 0;

  for (const auto& [section_name, dis] : target_discretizations_)
  {
    // Find out which section we have available for input. We can only do this on rank zero due
    // to large legacy sections that are not available everywhere. Communicate the result to all
    // ranks.
    std::map<std::string, bool> available_section;
    const int my_rank = Core::Communication::my_mpi_rank(comm_);
    if (my_rank == 0)
    {
      available_section[section_name + " ELEMENTS"] =
          input_.has_section(section_name + " ELEMENTS");
      available_section[section_name + " DOMAIN"] = input_.has_section(section_name + " DOMAIN");
      available_section[section_name + " GEOMETRY"] =
          input_.has_section(section_name + " GEOMETRY");
      Core::Communication::broadcast(available_section, 0, comm_);
    }
    else
    {
      Core::Communication::broadcast(available_section, 0, comm_);
    }

    const int num_sections_in_file =
        std::ranges::count_if(available_section, [](const auto& pair) { return pair.second; });
    if (num_sections_in_file > 1)
    {
      std::string found_sections;
      for (const auto& [section, exists] : available_section)
      {
        if (exists) found_sections += "'" + section + "' ";
      }
      FOUR_C_THROW(
          "Multiple options to read mesh for discretization '{}'. Only one is allowed.\n Found "
          "sections: {}",
          dis->name(), found_sections);
    }

    if (num_sections_in_file == 0)
    {
      // This used to be the default, so we use it for backwards compatibility.
      element_readers_.emplace_back(
          Core::IO::ElementReader(dis, input_, section_name + " ELEMENTS"));
      continue;
    }

    if (available_section[section_name + " ELEMENTS"])
    {
      element_readers_.emplace_back(
          Core::IO::ElementReader(dis, input_, section_name + " ELEMENTS"));
    }
    if (available_section[section_name + " DOMAIN"])
    {
      domain_readers_.emplace_back(Core::IO::DomainReader(dis, input_, section_name + " DOMAIN"));
    }
    if (available_section[section_name + " GEOMETRY"])
    {
      exodus_readers_.emplace_back(Internal::ExodusReader{
          .target_discretization = *dis,
          .section_name = section_name + " GEOMETRY",
      });
    }
  }

  graph_.resize(element_readers_.size());

  if (exodus_readers_.size() > 0)
  {
    read_mesh_from_exodus();
  }
  else
  {
    read_mesh_from_dat_file(max_node_id);
    rebalance();
    create_inline_mesh(max_node_id);
  }

  // last check if there are enough nodes
  {
    int local_max_node_id = max_node_id;
    Core::Communication::max_all(&local_max_node_id, &max_node_id, 1, comm_);

    if (max_node_id > 0 && max_node_id < Core::Communication::num_mpi_ranks(comm_))
      FOUR_C_THROW("Bad idea: Simulation with {} procs for problem with {} nodes",
          Core::Communication::num_mpi_ranks(comm_), max_node_id);
  }
}


MPI_Comm Core::IO::MeshReader::get_comm() const { return comm_; }


const Core::IO::Exodus::Mesh* Core::IO::MeshReader::get_exodus_mesh_on_rank_zero() const
{
  if (exodus_readers_.size() == 1)
  {
    return exodus_readers_.front().mesh_on_rank_zero.get();
  }
  return nullptr;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::read_mesh_from_dat_file(int& max_node_id)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::IO::MeshReader::read_mesh_from_dat_file");

  // read element information
  for (auto& element_reader : element_readers_) element_reader.read_and_distribute();

  // read nodes based on the element information
  read_nodes(input_, node_section_name_, element_readers_, max_node_id);
}


namespace
{
  std::pair<std::shared_ptr<Core::LinAlg::Map>, std::shared_ptr<Core::LinAlg::Map>>
  do_rebalance_discretization(const std::shared_ptr<const Core::LinAlg::Graph>& graph,
      Core::FE::Discretization& discret, Core::Rebalance::RebalanceType rebalanceMethod,
      Teuchos::ParameterList& rebalanceParams,
      const Core::IO::MeshReader::MeshReaderParameters& parameters, MPI_Comm comm)
  {
    std::shared_ptr<Core::LinAlg::Map> rowmap, colmap;

    switch (rebalanceMethod)
    {
      case Core::Rebalance::RebalanceType::hypergraph:
      {
        if (!Core::Communication::my_mpi_rank(comm))
          std::cout << "Redistributing using hypergraph .........\n";

        rebalanceParams.set("algorithm", "phg");
        
        std::tie(rowmap, colmap) = Core::Rebalance::rebalance_node_maps(*graph, rebalanceParams);
        break;
      }
      case Core::Rebalance::RebalanceType::multijagged:
      {
        if (!Core::Communication::my_mpi_rank(comm))
          std::cout << "Redistributing using multijagged .........\n";

        rebalanceParams.set("algorithm", "multijagged");

        rowmap = std::make_shared<Core::LinAlg::Map>(-1, graph->row_map().NumMyElements(),
            graph->row_map().MyGlobalElements(), 0, Core::Communication::as_epetra_comm(comm));
        colmap = std::make_shared<Core::LinAlg::Map>(-1, graph->col_map().NumMyElements(),
            graph->col_map().MyGlobalElements(), 0, Core::Communication::as_epetra_comm(comm));

        discret.redistribute(*rowmap, *colmap,
            {.assign_degrees_of_freedom = false,
                .init_elements = false,
                .do_boundary_conditions = false});

        std::shared_ptr<Core::LinAlg::MultiVector<double>> coordinates =
            discret.build_node_coordinates();

        std::tie(rowmap, colmap) = Core::Rebalance::rebalance_node_maps(
            *graph, rebalanceParams, nullptr, nullptr, coordinates);
        break;
      }
      case Core::Rebalance::RebalanceType::monolithic:
      {
        if (!Core::Communication::my_mpi_rank(comm))
          std::cout << "Redistributing using monolithic hypergraph .........\n";

        rebalanceParams.set("algorithm", "phg");

        rowmap = std::make_shared<Core::LinAlg::Map>(-1, graph->row_map().NumMyElements(),
            graph->row_map().MyGlobalElements(), 0, Core::Communication::as_epetra_comm(comm));
        colmap = std::make_shared<Core::LinAlg::Map>(-1, graph->col_map().NumMyElements(),
            graph->col_map().MyGlobalElements(), 0, Core::Communication::as_epetra_comm(comm));

        discret.redistribute(*rowmap, *colmap, {.do_boundary_conditions = false});

        std::shared_ptr<const Core::LinAlg::Graph> enriched_graph =
            Core::Rebalance::build_monolithic_node_graph(
                discret, Core::GeometricSearch::GeometricSearchParams(
                             parameters.geometric_search_parameters, parameters.io_parameters));

        std::tie(rowmap, colmap) =
            Core::Rebalance::rebalance_node_maps(*enriched_graph, rebalanceParams);
        break;
      }
      default:
        FOUR_C_THROW("Appropriate partitioning has to be set!");
    }

    return {rowmap, colmap};
  }

  void rebalance_discretization(Core::FE::Discretization& discret,
      const Core::LinAlg::Map& row_elements,
      const Core::IO::MeshReader::MeshReaderParameters& parameters, MPI_Comm comm)
  {
    std::shared_ptr<const Core::LinAlg::Graph> graph = nullptr;

    // Skip building the node graph if there are no elements
    if (row_elements.NumGlobalElements() > 0)
      graph = Core::Rebalance::build_graph(discret, row_elements);

    // Create partitioning parameters
    const double imbalance_tol =
        parameters.mesh_partitioning_parameters.get<double>("IMBALANCE_TOL");

    Teuchos::ParameterList rebalanceParams;
    rebalanceParams.set<std::string>("imbalance_tolerance", std::to_string(imbalance_tol));

    const int minele_per_proc =
        parameters.mesh_partitioning_parameters.get<int>("MIN_ELE_PER_PROC");
    const int max_global_procs = Core::Communication::num_mpi_ranks(comm);
    int min_global_procs = max_global_procs;

    if (minele_per_proc > 0) min_global_procs = row_elements.NumGlobalElements() / minele_per_proc;
    const int num_procs = std::min(max_global_procs, min_global_procs);
    rebalanceParams.set<std::string>("num_global_parts", std::to_string(num_procs));

    const auto rebalanceMethod = Teuchos::getIntegralValue<Core::Rebalance::RebalanceType>(
        parameters.mesh_partitioning_parameters, "METHOD");

    if (!Core::Communication::my_mpi_rank(comm))
      std::cout << "\nNumber of procs used for redistribution: " << num_procs << "\n";

    std::shared_ptr<Core::LinAlg::Map> rowmap, colmap;

    if (graph)
    {
      std::tie(rowmap, colmap) = do_rebalance_discretization(
          graph, discret, rebalanceMethod, rebalanceParams, parameters, comm);
    }
    else
    {
      rowmap = colmap = std::make_shared<Core::LinAlg::Map>(
          -1, 0, nullptr, 0, Core::Communication::as_epetra_comm(comm));
    }

    auto options_redistribution = Core::FE::OptionsRedistribution();
    if (rebalanceMethod == Core::Rebalance::RebalanceType::monolithic)
      options_redistribution.do_extended_ghosting = true;

    options_redistribution.assign_degrees_of_freedom = false;
    options_redistribution.init_elements = false;
    options_redistribution.do_boundary_conditions = false;

    discret.redistribute(*rowmap, *colmap, options_redistribution);

    Core::Rebalance::Utils::print_parallel_distribution(discret);
  }
}  // namespace

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::rebalance() const
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::IO::MeshReader::Rebalance");

  for (const auto& element_reader : element_readers_)
  {
    rebalance_discretization(
        *element_reader.get_dis(), *element_reader.get_row_elements(), parameters_, comm_);
  }
}


void Core::IO::MeshReader::read_mesh_from_exodus()
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::IO::MeshReader::read_mesh_from_exodus");

  auto my_rank = Core::Communication::my_mpi_rank(comm_);

  for (auto& exodus_reader : exodus_readers_)
  {
    // We cannot create the map right away. First, we need to figure out how many elements there
    // are. Since the code is rather different on rank 0 and other ranks, we will set this pointer
    // to nullptr and create it later.
    std::unique_ptr<Core::LinAlg::Map> linear_element_map;

    // All the work is done on rank 0. The other ranks will receive the data.
    if (my_rank == 0)
    {
      if (!input_.has_section(exodus_reader.section_name)) return;

      InputParameterContainer data;
      input_.match_section(exodus_reader.section_name, data);

      const auto& geometry_data = data.group(exodus_reader.section_name);
      const auto& exodus_file = geometry_data.get<std::filesystem::path>("FILE");
      exodus_reader.mesh_on_rank_zero = std::make_unique<Exodus::Mesh>(exodus_file.string());
      const auto& mesh = *exodus_reader.mesh_on_rank_zero;

      // Initial implementation:
      // - read all information on rank 0, construct discretization, rebalance afterwards

      const auto& element_block_data = geometry_data.get_list("ELEMENT_BLOCKS");

      // Exodus does not number elements internally. We number the elements in the order we
      // encounter them in the blocks.
      int ele_count = 0;
      for (const auto& [eb_id, eb] : mesh.get_element_blocks())
      {
        // Look into the input file to find out which elements we need to assign to this block.
        const int eb_id_copy = eb_id;  // work around compiler warning in clang18
        auto current_block_data = std::ranges::find_if(element_block_data,
            [eb_id_copy](const auto& e) { return e.template get<int>("ID") == eb_id_copy; });
        FOUR_C_ASSERT_ALWAYS(current_block_data != element_block_data.end(),
            "The exodus file contains element block {} but you did not specify information for "
            "this block in the input file.",
            eb_id);

        const auto& element_name = current_block_data->get<std::string>("ELEMENT_NAME");
        const auto cell_type = Exodus::shape_to_cell_type(eb.get_shape());
        const auto cell_type_string = FE::cell_type_to_string(cell_type);

        Core::Elements::ElementDefinition ed;
        ed.setup_valid_element_lines();
        const auto& linedef = ed.element_lines(element_name, cell_type_string);

        // The spec for elements also contains the nodes for the legacy input.
        // Thus, we fake a string here that contains the cell_type followed by the appropriate
        // number of dummy nodes (that we are not going to use).
        std::stringstream ss;
        ss << cell_type_string;
        const int numnodes = FE::cell_type_switch(
            cell_type, [](auto cell_type_t) { return FE::num_nodes<cell_type_t()>; });
        for (int i = 0; i < numnodes; ++i) ss << " " << 0;  // dummy node id
        ss << " " << current_block_data->get<std::string>("ELEMENT_DATA");
        std::string element_string = ss.str();

        Core::IO::ValueParser element_parser{
            element_string, {.user_scope_message = "While reading element data: "}};
        Core::IO::InputParameterContainer element_data;
        linedef.fully_parse(element_parser, element_data);


        for (const auto& ele_nodes : *eb.get_ele_conn() | std::views::values)
        {
          auto ele = Core::Communication::factory(element_name, cell_type_string, ele_count, 0);
          if (!ele) FOUR_C_THROW("element creation failed");
          ele->set_node_ids(ele_nodes.size(), ele_nodes.data());
          ele->read_element(element_name, cell_type_string, element_data);
          exodus_reader.target_discretization.add_element(ele);

          ele_count++;
        }
      }

      int num_ele = mesh.get_num_ele();
      Core::Communication::broadcast(num_ele, 0, comm_);
      linear_element_map = std::make_unique<Core::LinAlg::Map>(
          num_ele, 0, Core::Communication::as_epetra_comm(comm_));

      std::vector<int> gid_list(mesh.get_num_ele());
      std::iota(gid_list.begin(), gid_list.end(), 0);
      exodus_reader.target_discretization.proc_zero_distribute_elements_to_all(
          *linear_element_map, gid_list);

      // Now add all the nodes to the discretization on rank 0. They are distributed later during
      // the rebalancing process.
      for (const auto& [id, coords] : mesh.get_nodes())
      {
        auto node = std::make_shared<Core::Nodes::Node>(id, coords, 0);
        exodus_reader.target_discretization.add_node(node);
      }
    }
    // Other ranks
    else
    {
      int num_ele;
      Core::Communication::broadcast(num_ele, 0, comm_);
      linear_element_map = std::make_unique<Core::LinAlg::Map>(
          num_ele, 0, Core::Communication::as_epetra_comm(comm_));

      std::vector<int> gid_list;
      exodus_reader.target_discretization.proc_zero_distribute_elements_to_all(
          *linear_element_map, gid_list);
    }

    FOUR_C_ASSERT(linear_element_map, "Internal error: nullptr.");
    rebalance_discretization(
        exodus_reader.target_discretization, *linear_element_map, parameters_, comm_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::create_inline_mesh(int& max_node_id) const
{
  for (const auto& domain_reader : domain_readers_)
  {
    // communicate node offset to all procs
    int local_max_node_id = max_node_id;
    Core::Communication::max_all(&local_max_node_id, &max_node_id, 1, comm_);

    domain_reader.create_partitioned_mesh(max_node_id);
    domain_reader.complete();
    max_node_id = domain_reader.my_dis()->node_row_map()->MaxAllGID() + 1;
  }
}

FOUR_C_NAMESPACE_CLOSE
