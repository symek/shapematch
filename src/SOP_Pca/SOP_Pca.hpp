//
// Created by symek on 9/13/19.
//
#pragma once
#include <SOP/SOP_Node.h>

namespace pca {
    using Vertices    = Eigen::MatrixXd;
    void copy_position_to_eigen(const GU_Detail * gdp, Vertices & matrix)
    {
        // Vertices matrix;
        matrix.conservativeResize(3, gdp->getNumPoints());

        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {
            const UT_Vector3 pos = gdp->getPos3(ptoff);
            const GA_Index   idx = gdp->pointIndex(ptoff);
            matrix(0, idx) = pos.x();
            matrix(1, idx) = pos.y();
            matrix(2, idx) = pos.z();
        }
    }

    class SOP_Pca : public SOP_Node
    {
    public:
        SOP_Pca(OP_Network *net, const char *name, OP_Operator *op);
        virtual ~SOP_Pca();
//        virtual bool updateParmsFlags() override;
        static PRM_Template      myTemplateList[];
        static OP_Node      *myConstructor(OP_Network*, const char *,
                                           OP_Operator *);
    protected:
        /// Method to cook geometry for the SOP
        virtual OP_ERROR         cookMySop(OP_Context &context) override;

    };

} // End  namespace
