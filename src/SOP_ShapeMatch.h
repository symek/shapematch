//
// Created by symek on 9/12/19.
//


#pragma once
#include <SOP/SOP_Node.h>

namespace shapematch {

    typedef double Scalar;
    typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic, Eigen::ColMajor> Vertices;

    enum ALIGN_METHOD {
        RIGID,
        SPARSE_ICP,
        REWEIGHTED_ICP,
        INTEL_FGR,
        RIGID_CPD,
        NONRIGID_CPD,
    };

    enum WEIGHT_FUNC {
        PNORM,
        TUKEY,
        FAIR,
        LOGISTIC,
        TRIMMED,
        NONE,
    };

    bool copy_float_to_eigen(const GA_Attribute * attr, Eigen::VectorXd & weightsV)
    {
        const GA_Detail & gdp = attr->getDetail();
        GA_ROHandleF  weights_h(attr);

        if (weights_h.isValid())
        {
            weightsV.conservativeResize(gdp.getNumPoints());
            UT_ASSERT(weightsV.rows() == gdp.getNumPoints());
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(&gdp, ptoff)
            {
                const float w = weights_h.get(ptoff);
                const GA_Index idx = gdp.pointIndex(ptoff);
                weightsV(idx) = static_cast<double>(w);
            }
            return true;
        }
        return false;
    }

    void copy_position_to_eigen(const GU_Detail * gdp, Vertices & matrix)
    {
        // Vertices matrix;
        matrix.conservativeResize(Eigen::NoChange, gdp->getNumPoints());

        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {
            const UT_Vector3 pos = gdp->getPos3(ptoff);
            const GA_Index   idx = gdp->pointIndex(ptoff);
            matrix(0, idx) = pos.x();
            matrix(1, idx) = pos.y();
            matrix(2, idx) = pos.z();
        }
        // return matrix;
    }

    Eigen::MatrixXd copy_position_to_eigen_rows(const GU_Detail * gdp)
    {
        Eigen::MatrixXd matrix;
        matrix.resize(gdp->getNumPoints(), 3);

        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(gdp, ptoff)
        {
            const UT_Vector3 pos = gdp->getPos3(ptoff);
            const GA_Index   idx = gdp->pointIndex(ptoff);
            matrix(idx, 0) = static_cast<double>(pos.x());
            matrix(idx, 1) = static_cast<double>(pos.y());
            matrix(idx, 2) = static_cast<double>(pos.z());
        }
        return matrix;
    }

    class SOP_ShapeMatch : public SOP_Node
    {
    public:
        SOP_ShapeMatch(OP_Network *net, const char *name, OP_Operator *op);
        virtual ~SOP_ShapeMatch();
//        virtual bool updateParmsFlags() override;
        static PRM_Template      myTemplateList[];
        static OP_Node      *myConstructor(OP_Network*, const char *,
                                           OP_Operator *);
    protected:
        /// Method to cook geometry for the SOP
        virtual OP_ERROR         cookMySop(OP_Context &context) override;
    private:

        int     MAXITER(fpreal t)       { return evalInt("maxiterations", 0, t); }
        fpreal  STOPCRITERIA(fpreal t)  { return evalFloat("stopcritera", 0, t); }

        void    add_detail_array(const Eigen::MatrixXd & , const char* attr_name="rigid_xform");

    };

} // End pcallign namespace
