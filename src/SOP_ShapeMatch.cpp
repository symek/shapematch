//
// Created by symek on 9/12/19.
//

#include <UT/UT_DSOVersion.h>
#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_AutoLockInputs.h>
#include <SYS/SYS_Math.h>
#include <Eigen/Dense>

#include <cpd/nonrigid.hpp>
#include <cpd/gauss_transform.hpp>
#include <cpd/gauss_transform_fgt.hpp>

#include "SOP_ShapeMatch.h"

using namespace shapematch;


void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
            "shapematch",
            "ShapeMatch",
            SOP_ShapeMatch::myConstructor,
            SOP_ShapeMatch::myTemplateList,
            2,
            2,
            0));
}

static PRM_Name names[] = {
        PRM_Name("maxiterations",  "Max inner iteration"),
        PRM_Name("stopcritera",   "Stopping criteria"),
        PRM_Name("docorrespondance", "Compute correspondence"),
        PRM_Name("outliers",         "Outliers factor"),

};



static PRM_Default maxIterDefault(10);
static PRM_Default stopDefault(0.3);
static PRM_Default outliersDefault(0.5);

PRM_Template
        SOP_ShapeMatch::myTemplateList[] = {
        PRM_Template(PRM_INT_J, 1, &names[0], &maxIterDefault), // max iterations (rigid, sparse, reweighted sparse)
        PRM_Template(PRM_FLT_LOG, 1, &names[1], &stopDefault),    // stop criteria
        PRM_Template(PRM_TOGGLE,1, &names[2], PRMzeroDefaults), // compute correspondance rigid CPD
        PRM_Template(PRM_FLT_LOG, 1,&names[3], &outliersDefault),    // outliers CPD
        PRM_Template(),
};




OP_Node *
SOP_ShapeMatch::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_ShapeMatch(net, name, op);
}

SOP_ShapeMatch::SOP_ShapeMatch(OP_Network *net, const char *name, OP_Operator *op)
        : SOP_Node(net, name, op)
{
    mySopFlags.setManagesDataIDs(true);
}

SOP_ShapeMatch::~SOP_ShapeMatch() {}

OP_ERROR
SOP_ShapeMatch::cookMySop(OP_Context &context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    duplicatePointSource(0, context);

    // Get second geometry:
    const GU_Detail * target_gdp = inputGeo(1);

    // any points?:
    if (target_gdp->getNumPoints() == 0 || gdp->getNumPoints() == 0) {
        addError(SOP_MESSAGE, "Needs two points clouds to align.");
        return error();
    }


    const int   max_iterations       =  evalInt("maxiterations", 0, t);
    const float stop_critera         =  evalFloat("stopcritera", 0, t);
    const int   allow_correspondance = evalInt("docorrespondance", 0, t);
    const float outliers             = evalFloat("outliers", 0, t);

    Vertices source, target;
    copy_position_to_eigen(gdp, source);
    copy_position_to_eigen(target_gdp, target);

    UT_ASSERT(gdp->getNumPoints() == source.cols());
    UT_ASSERT(target_gdp->getNumPoints() == target.cols());
    UT_ASSERT(source.cols() == target.cols());


    Eigen::MatrixXd sourcet = source.transpose();
    Eigen::MatrixXd targett = target.transpose();

    cpd::Nonrigid nonrigid;

    nonrigid.correspondence(static_cast<bool>(allow_correspondance));
    nonrigid.outliers(outliers);
    nonrigid.max_iterations(SYSmax(max_iterations,1));
    nonrigid.tolerance(stop_critera);

    // nonrigid.gauss_transform(std::move(
    // std::unique_ptr<cpd::GaussTransform>(new cpd::GaussTransformFgt())));

    cpd::NonrigidResult result = nonrigid.run(sourcet, targett);

    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        const GA_Index i = gdp->pointIndex(ptoff);
        const UT_Vector3 pos(result.points(i, 0), result.points(i, 1), result.points(i, 2));
        gdp->setPos3(ptoff, pos);
    }

    gdp->getP()->bumpDataId();
    return error();
}