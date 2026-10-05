#pragma once

#include "ortools/sat/cp_model.h"
#include <vector>
#include "safe_vector.h"

using ORIntervalVars = Vector<operations_research::sat::IntervalVar>;
using ORIntVars1D = Vector<operations_research::sat::IntVar>;
using ORIntVars2D = Vector<Vector<operations_research::sat::IntVar>>;
using ORBoolVars1D = Vector<operations_research::sat::BoolVar>;
using ORBoolVars2D = Vector<Vector<operations_research::sat::BoolVar>>;
using ORBoolVars3D = Vector<Vector<Vector<operations_research::sat::BoolVar>>>;
using ORLinExpr = operations_research::sat::LinearExpr;
