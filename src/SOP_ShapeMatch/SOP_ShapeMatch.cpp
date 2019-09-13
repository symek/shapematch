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
#include <cpd/gauss_transform_fgt.hpp>

#include "SOP_ShapeMatch.h"

using namespace shapematch;


void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
            "shapematch",
            "Shape Match",
            SOP_ShapeMatch::myConstructor,
            SOP_ShapeMatch::myTemplateList,
            2,
            2,
            0));
}

static PRM_Name names[] = {
        PRM_Name("maxiterations",    "Max iteration"),
        PRM_Name("tolerance",        "Tolerance "),
        PRM_Name("outliers",         "Outliers factor"),
        PRM_Name("usefgt",            "Use FastGaussTransform"),
        PRM_Name("fgtmethod",         "FGT Method"),
        PRM_Name("breakpoint",        "Breakpoint"),
        PRM_Name("docorrespondence", "Compute correspondence"),
        PRM_Name("applyshapematch",   "Apply Shape Match"),

};

static PRM_Name  fgtMethodChoices[] =
        {
                PRM_Name("0", "Direct Tree"),
                PRM_Name("1", "Improved FGT"),
                PRM_Name("2", "Switched"),
                PRM_Name(0)
        };

static PRM_ChoiceList  methodMenu(PRM_CHOICELIST_SINGLE, fgtMethodChoices);

static PRM_Default maxIterDefault(16);
static PRM_Default toleranceDefault(1e-5);
static PRM_Default outliersDefault(1.0);
static PRM_Default breakpointDefault(0.2);

PRM_Template
        SOP_ShapeMatch::myTemplateList[] = {
        PRM_Template(PRM_INT_J, 1,   &names[0],    &maxIterDefault), // max iterations
        PRM_Template(PRM_FLT_LOG, 1, &names[1],    &toleranceDefault),    // stop criteria
        PRM_Template(PRM_FLT_LOG, 1, &names[2],    &outliersDefault),    // outliers
        PRM_Template(PRM_TOGGLE,1,   &names[3],    PRMoneDefaults), // use fgt
        PRM_Template(PRM_ORD,   1,   &names[4], 0, &methodMenu), // fgt method
        PRM_Template(PRM_FLT_LOG, 1, &names[5],    &breakpointDefault),    // breakpoint
        PRM_Template(PRM_TOGGLE,1,   &names[6],    PRMzeroDefaults), // compute correspondance
//        PRM_Template(PRM_FLT_LOG, 1, &names[7],    PRMoneDefaults),    // breakpoint
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

    // any points?:
    if (target_gdp->getNumPoints() != gdp->getNumPoints()) {
        addError(SOP_MESSAGE, "Point counts should match (for now).");
        return error();
    }


    const int   max_iterations       = evalInt("maxiterations", 0, t);
    const float tolerance            = evalFloat("tolerance", 0, t);
    const float outliers             = evalFloat("outliers", 0, t);
    const int   use_fgt              = evalInt("usefgt", 0, t);
    const int   fgt_method           = evalInt("fgtmethod", 0, t);
    const float breakpoint           = evalFloat("breakpoint", 0, t);
    const int   do_correspondence    = evalInt("docorrespondence", 0, t);

    Vertices source, target;
    copy_position_to_eigen(gdp, source);
    copy_position_to_eigen(target_gdp, target);

    UT_ASSERT(gdp->getNumPoints() == source.cols());
    UT_ASSERT(target_gdp->getNumPoints() == target.cols());
    UT_ASSERT(source.cols() == target.cols());

    Eigen::MatrixXd sourcet = source.transpose();
    Eigen::MatrixXd targett = target.transpose();

    cpd::Nonrigid nonrigid;

    nonrigid.max_iterations(SYSmax(max_iterations,1));
    nonrigid.outliers(outliers);
    nonrigid.correspondence(static_cast<bool>(do_correspondence));
    nonrigid.tolerance(tolerance);

    if (use_fgt) {
        std::unique_ptr<cpd::GaussTransformFgt> fgt(new cpd::GaussTransformFgt());
        fgt->method(static_cast<cpd::FgtMethod>(fgt_method));
        fgt->breakpoint(breakpoint);
        nonrigid.gauss_transform(std::move(fgt));
    }

    cpd::NonrigidResult result = nonrigid.run(sourcet, targett);
    auto runtime = result.runtime;
    //std::chrono::duration_cast<std::chrono::seconds>(result.runtime);
    std::string message("Compute time: ");
    message += std::to_string(runtime.count()/1E6);
    message += "s, iterations: " + std::to_string(result.iterations);
    addMessage(SOP_MESSAGE, message.c_str());

    GA_Attribute * correspondence_attr = nullptr;
    if (do_correspondence) {
        correspondence_attr = gdp->addFloatTuple(GA_ATTRIB_POINT, "pointmatch", 1);
        if (!correspondence_attr) {
            addError(SOP_MESSAGE, "Can't create 'pointmatch' attribute.");
            return error();
        }
    }
    GA_RWHandleI correspondence_h(correspondence_attr);
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        const GA_Index i = gdp->pointIndex(ptoff);
        const UT_Vector3 pos(result.points(i, 0), result.points(i, 1), result.points(i, 2));
        gdp->setPos3(ptoff, pos);

        if (do_correspondence) {
            correspondence_h.set(ptoff, result.correspondence(i));
        }
    }

    gdp->getP()->bumpDataId();
    return error();
}