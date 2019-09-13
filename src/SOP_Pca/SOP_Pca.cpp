//
// Created by symek on 9/13/19.
//

#include <UT/UT_DSOVersion.h>
#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_AutoLockInputs.h>

#include <Eigen/Dense>

#include "SOP_Pca.hpp"
#include "pca.hpp"

using namespace pca;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
            "principlecomponent",
            "Priniciple Component",
            SOP_Pca::myConstructor,
            SOP_Pca::myTemplateList,
            1,
            1,
            0));
}

static PRM_Name names[] = {
        PRM_Name("variance",        "Variance "),
        PRM_Name("shift",         "Shift Data"),
        PRM_Name("ortho",            "Orthogonalize PCA"),


};


static PRM_Default varianceDefault(0.01);


PRM_Template
        SOP_Pca::myTemplateList[] = {
        PRM_Template(PRM_FLT_LOG, 1, &names[0],    &varianceDefault),    // stop criteria
        PRM_Template(PRM_TOGGLE,1,   &names[1],    PRMzeroDefaults), // use fgt
        PRM_Template(PRM_TOGGLE,1,   &names[2],    PRMzeroDefaults), // use fgt

        PRM_Template(),
};

OP_Node *
SOP_Pca::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_Pca(net, name, op);
}

SOP_Pca::SOP_Pca(OP_Network *net, const char *name, OP_Operator *op)
        : SOP_Node(net, name, op)
{
    mySopFlags.setManagesDataIDs(true);
}

SOP_Pca::~SOP_Pca() {}

OP_ERROR
SOP_Pca::cookMySop(OP_Context &context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    fpreal t = context.getTime();
    duplicatePointSource(0, context);

    // any points?:
    if (gdp->getNumPoints() == 0) {
        addError(SOP_MESSAGE, "Needs some points to work with.");
        return error();
    }

//
    const float variance      = evalFloat("variance", 0, t);
    const int   shift         = evalInt("shift", 0, t);
    const int   orthogonalize = evalInt("ortho", 0, t);


    Vertices source;
    copy_position_to_eigen(gdp, source);

    UT_ASSERT(gdp->getNumPoints() == source.cols());

    Vertices pcamatrix;
    computePCA(source, pcamatrix, variance, static_cast<bool>(shift), static_cast<bool>(orthogonalize));

    const GA_Range point_range = gdp->getPointRangeSlice(0, -1);
    gdp->destroyPointOffsets(point_range, GA_Detail::GA_DESTROY_DEGENERATE);
//
    gdp->appendPointBlock(pcamatrix.cols());
//
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        const GA_Index index = gdp->pointIndex(ptoff);
        UT_Vector3 pos(pcamatrix(0, index), pcamatrix(1, index), pcamatrix(2, index));
        gdp->setPos3(ptoff, pos);
    }

    gdp->getP()->bumpDataId();
    return error();
}