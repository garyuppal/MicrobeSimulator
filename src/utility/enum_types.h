#ifndef MICROBESIMULATOR_ENUM_TYPES_H
#define MICROBESIMULATOR_ENUM_TYPES_H

#include <string>

// boost library:
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

namespace MicrobeSimulator{

	enum class RunMode : int
	{
		FEM_CG = 0, FEM_DG = 1, FDM = 2 //, DEBUGGING = 3
	};

	enum class SimulatorType : int
	{
		BACTERIA, INTCHT, AGING	, TWOPG
	};

	// enum class BoundaryCondition : int
	// {
	// 	WRAP = 0, REFLECT = 1, OPEN = 2 // open end will be upper end for now
	// }; // only makes sense for hyper cube
	
  	enum class VelocityType : int
	{
		NO_FLOW = 0, STOKES = 1, CONSTANT = 2, 
		COUETTE = 3, SQUARE_PIPE = 4, VORTEX = 5
	}; // velocity type


	// @todo : maybe change name to geometry type ...

	enum class SourceImplementation : int
	{
		NATIVE = 0, POINT_CELL = 1, SOURCE_MAP = 2
	};

	// enum class SimulationType : int
	// {
	// 	PIPE = 1, MIXER = 2, FILTER = 3
	// }; // should know from mesh type

	SimulatorType stringToSimulatorType(std::string value)
	{
		if( boost::iequals(value, "Bacteria") )
			return SimulatorType::BACTERIA;
		else if( boost::iequals(value, "Ic_as") ||  boost::iequals(value, "INTCHT"))
			return SimulatorType::INTCHT;
		else if( boost::iequals(value, "Aging"))
			return SimulatorType::AGING;

		// legacy:
		else if( boost::iequals(value,"TWOPG") )
			return SimulatorType::TWOPG;
		
		else
		{
			std::string msg = "Invalid string! Could not convert <"
							+ value
							+ "> to run mode type.";
			throw std::invalid_argument(msg.c_str());
		}
			// throw std::invalid_argument("Invalid string! Could not convert to simulator type");
	}

	std::string getSimulatorTypeString(SimulatorType sim_type)
	{
		std::string result;
		if(sim_type == SimulatorType::BACTERIA)
			result = "Bacteria Simulation";
		else if(sim_type == SimulatorType::INTCHT)
			result = "Intermittent Cheating Simulation";
		else if(sim_type == SimulatorType::AGING)
			result = "Aging simulation";
		else if(sim_type == SimulatorType::TWOPG)
			result = "Two public good simulation";
		else 
			result = "ERROR";

		return result;
	}

	// run mode -> convert to implementation method (each chemical can 
	// be implemented differently as needed )
	std::string getRunModeString(RunMode run_mode)
	{
		std::string result = "";
		if(run_mode == RunMode::FEM_CG)
			result = "Finite element (CG)";
		else if(run_mode == RunMode::FEM_DG)
			result = "Discontinuous Galerkin finite element (DG)";
		else if(run_mode == RunMode::FDM)
			result = "Finite difference - Forward Euler";
		else 
			result = "ERROR";

		return result;
	}


	RunMode stringToRunMode(std::string value)
	{
		if( boost::iequals(value,"FEM_CG") )
			return RunMode::FEM_CG;
		else if( boost::iequals(value,"FEM_DG") )
			return RunMode::FEM_DG;
		else if( boost::iequals(value,"FDM") )
			return RunMode::FDM;
		else
		{
			std::string msg = "Invalid string! Could not convert <"
							+ value
							+ "> to run mode type.";
			throw std::invalid_argument(msg.c_str());
		}
	}	


	std::string getVelocityTypeString(VelocityType vtype)
	{
		std::string result = "";
		if(vtype == VelocityType::NO_FLOW)
			result = "No flow";
		else if(vtype == VelocityType::STOKES)
			result = "Numerical (stokes) flow field";
		else if(vtype == VelocityType::CONSTANT)
			result = "Constant flow";
		else if(vtype == VelocityType::COUETTE)
			result = "Planar Couette flow";
		else if(vtype == VelocityType::SQUARE_PIPE)
			result = "Hagen-Poiseuille square pipe flow";
		else if(vtype == VelocityType::VORTEX)
			result = "Rankine vortex flow";
		else
			result = "ERROR";

		return result;
	}


	VelocityType stringToVelocityType(std::string value)
	{												// legacy v
		if( (boost::iequals(value, "NO FLOW")) || (boost::iequals(value, "NO FLOW")) )
			return VelocityType::NO_FLOW;
		else if( boost::iequals(value, "STOKES") )
			return VelocityType::STOKES;
		else if( boost::iequals(value, "CONSTANT") )
			return VelocityType::CONSTANT;
		else if( boost::iequals(value, "COUETTE") )
			return VelocityType::COUETTE;
		else if( boost::iequals(value, "SQUARE PIPE") )
			return VelocityType::SQUARE_PIPE;
		else if( boost::iequals(value, "VORTEX") )
			return VelocityType::VORTEX;
		else
		{
			std::string msg = "Invalid string! Could not convert <"
							+ value
							+ "> to velocity type.";
			throw std::invalid_argument(msg.c_str());
		}
			// throw std::invalid_argument("Invalid string! Could not convert to velocity type");
	}




	// enum class MeshType : int
	// {
	// 	FILE_MESH = 0, BOX_MESH = 1, SQUARE_CHEESE = 2, HEX_CHEESE = 3,
	// 	MIXER = 4
	// };



	std::string getSourceImplementationString(SourceImplementation s_imp)
	{
		std::string result = "";
		if(s_imp == SourceImplementation::NATIVE)
			result = "Native";
		else if(s_imp == SourceImplementation::POINT_CELL)
			result = "Point cell field";
		else if(s_imp == SourceImplementation::SOURCE_MAP)
			result = "Source map";
		else
			result = "ERROR";

		return result;
	}


	SourceImplementation stringToSourceImplementation(std::string value)
	{
		if( boost::iequals(value,"NATIVE") )
			return SourceImplementation::NATIVE;
		else if( boost::iequals(value,"POINT_CELL") )
			return SourceImplementation::POINT_CELL;
		else if( boost::iequals(value,"SOURCE_MAP") )
			return SourceImplementation::SOURCE_MAP;
		else
			throw std::invalid_argument("Invalid string! Could not convert to source implementation type");	
	}
}
#endif

