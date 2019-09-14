//
// Created by symek on 9/14/19.
//

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include "converters.hpp"


namespace SOP_IGL {


    void detail_to_eigen(const GU_Detail &gdp, Eigen::MatrixXd &points, Eigen::MatrixXi &faces)
    {
        points.resize(0,0); faces.resize(0,0);
        points.resize(gdp.getNumPoints(), 3);
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(&gdp, ptoff) {
                const UT_Vector3 pos = gdp.getPos3(ptoff);
                points(static_cast<uint>(ptoff), 0) = pos.x();
                points(static_cast<uint>(ptoff), 1) = pos.y();
                points(static_cast<uint>(ptoff), 2) = pos.z();
            }

        UT_Array<const GA_Primitive *> prims;
        const GA_PrimitiveTypeId type(GA_PRIMPOLY);
        gdp.getPrimitivesOfType(type, prims);
        faces.resize(prims.size(), 3);
        UT_Array<const GA_Primitive*>::const_iterator it;
        for (it = prims.begin(); !it.atEnd(); ++it)
        {
            const GEO_Primitive *prim = static_cast<const GEO_Primitive*>(*it);
            const int vertex_count    = prim->getVertexCount();
            for (int vt=0; vt<vertex_count; ++vt) {
                const GA_Offset voff = prim->getPointOffset(vt);
                faces(prim->getMapIndex(), SYSmin(vt, 2)) = static_cast<int>(voff);
            }
        }
    }

    void detail_to_eigen(const GU_Detail &gdp, Eigen::MatrixXd &points,\
                     Eigen::MatrixXi &faces, Eigen::MatrixXi &tetras)
    {
        points.resize(0,0); faces.resize(0,0); tetras.resize(0,0);
        points.resize(gdp.getNumPoints(), 3);
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(&gdp, ptoff) {
                const UT_Vector3 pos = gdp.getPos3(ptoff);
                points(static_cast<uint>(ptoff), 0) = pos.x();
                points(static_cast<uint>(ptoff), 1) = pos.y();
                points(static_cast<uint>(ptoff), 2) = pos.z();
            }

        {
            UT_Array<const GA_Primitive *> prims;
            const GA_PrimitiveTypeId type(GA_PRIMPOLY);
            gdp.getPrimitivesOfType(type, prims);
            faces.resize(prims.size(), 3);
            UT_Array<const GA_Primitive*>::const_iterator it;
            for (it = prims.begin(); !it.atEnd(); ++it) {
                const GEO_Primitive *prim = static_cast<const GEO_Primitive*>(*it);
                int primitive_index = 0;
                const int vertex_count = prim->getVertexCount();
                for (int vt=0; vt<vertex_count; ++vt) {
                    const GA_Offset voff = prim->getPointOffset(vt);
                    faces(primitive_index, SYSmin(vt, 2)) = static_cast<int>(voff);
                    primitive_index++;
                }
            }
        }

        {
            UT_Array<const GA_Primitive *> prims;
            const GA_PrimitiveTypeId type(GA_PRIMTETRAHEDRON);
            gdp.getPrimitivesOfType(type, prims);
            tetras.resize(prims.size(), 4);
            UT_Array<const GA_Primitive*>::const_iterator it;
            for (it = prims.begin(); !it.atEnd(); ++it) {
                const GEO_Primitive *prim = static_cast<const GEO_Primitive*>(*it);
                int primitive_index = 0;
                const int vertex_count = prim->getVertexCount();
                for (int vt=0; vt<vertex_count; ++vt) {
                    const GA_Offset voff = prim->getPointOffset(vt);
                    tetras(primitive_index, SYSmin(vt, 3)) = static_cast<int>(voff);
                    primitive_index++;
                }
            }
        }
    }

    void eigen_to_detail(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, GU_Detail &gdp)
    {

        GA_Offset ptoff;
        UT_Vector3 pos;
        gdp.appendPointBlock((GA_Size)points.rows());
        GA_FOR_ALL_PTOFF(&gdp, ptoff) {
                pos.x() = points(static_cast<uint>(ptoff), 0);
                pos.y() = points(static_cast<uint>(ptoff), 1);
                pos.z() = points(static_cast<uint>(ptoff), 2);
                gdp.setPos3(ptoff, pos);
            }

        // TODO: optimize!
        for (uint i = 0; i < faces.rows(); ++i) {
            GU_PrimPoly *prim = GU_PrimPoly::build(&gdp, 0, false, false);
            prim->appendVertex((GA_Index)faces(i, 0));
            prim->appendVertex((GA_Index)faces(i, 1));
            prim->appendVertex((GA_Index)faces(i, 2));
        }
    }

    void eigen_to_detail_points(const Eigen::MatrixXd &points, GU_Detail &gdp)
    {
        GA_Offset ptoff;
        GA_FOR_ALL_PTOFF(&gdp, ptoff)
            {
                GA_Index ptidx = gdp.pointIndex(ptoff);
                if ((uint)ptidx < points.rows()) {
                    UT_Vector3 pos(points((uint)ptidx, 0),
                                   points((uint)ptidx, 1),
                                   points((uint)ptidx, 2));
                    gdp.setPos3(ptoff, pos);
                }

            }
    }

    void eigen_to_point_attribF(const Eigen::MatrixXd &M, GU_Detail *gdp, const char* attribName)
    {
        const uint attrib_size = M.cols();
        if (attrib_size != 1)
            return; // FIXME: silent quit not very elegant

        GA_RWHandleF  attrib_h(gdp->addFloatTuple(GA_ATTRIB_POINT, attribName, 1));
        if (attrib_h.isValid()) {
            GA_Offset ptoff;
            // UT_ASSERT(M.rows() == gdp->getNumPoints());
            GA_FOR_ALL_PTOFF(gdp, ptoff) {
                    GA_Index ptidx = gdp->pointIndex(ptoff);
                    if ((uint)ptidx < M.rows()) {
                        const float val = M((uint)ptidx, 0);
                        attrib_h.set(ptoff, val);
                    }
                }
        }

    }

    void eigen_to_point_attribV(const Eigen::MatrixXd &M, GU_Detail *gdp, const char* attribName)
    {
        const uint attrib_size = M.cols();
        if (attrib_size != 3)
            return;

        GA_RWHandleV3  attrib_h(gdp->addFloatTuple(GA_ATTRIB_POINT, attribName, 3));
        if (attrib_h.isValid()) {
            GA_Offset ptoff;
            // UT_ASSERT(M.rows() == gdp->getNumPoints());
            GA_FOR_ALL_PTOFF(gdp, ptoff) {
                    GA_Index ptidx = gdp->pointIndex(ptoff);
                    if ((uint)ptidx < M.rows()) {
                        UT_Vector3 val(M((uint)ptidx, 0),
                                       M((uint)ptidx, 1),
                                       M((uint)ptidx, 2));
                        attrib_h.set(ptoff, val);
                    }
                }
        }

    }

} // end of namespace SOP_IGL
