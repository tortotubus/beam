#include "DLMInternal.H"

namespace ELFF::Fluid::DLM::Detail {

std::array<amrex::Real, 3>&
markerField(MarkerPoint& marker, MarkerField field)
{
    switch (field) {
    case MarkerField::Lambda:
        return marker.lambda;
    case MarkerField::Residual:
        return marker.residual;
    case MarkerField::Search:
        return marker.search;
    case MarkerField::Response:
        return marker.response;
    }
    return marker.residual;
}

std::array<amrex::Real, 3> const&
markerField(MarkerPoint const& marker, MarkerField field)
{
    switch (field) {
    case MarkerField::Lambda:
        return marker.lambda;
    case MarkerField::Residual:
        return marker.residual;
    case MarkerField::Search:
        return marker.search;
    case MarkerField::Response:
        return marker.response;
    }
    return marker.residual;
}

bool
hasNonzeroMarkerField(std::vector<MarkerPoint> const& markers, MarkerField field)
{
    for (MarkerPoint const& marker : markers) {
        auto const& values = markerField(marker, field);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (std::abs(values[dir]) > 0.0) {
                return true;
            }
        }
    }
    return false;
}

void
copyMarkerField(std::vector<MarkerPoint>& markers, MarkerField dst, MarkerField src)
{
    for (MarkerPoint& marker : markers) {
        auto& dst_values = markerField(marker, dst);
        auto const& src_values = markerField(marker, src);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            dst_values[dir] = src_values[dir];
        }
    }
}

void
resetMarkerUzawaFields(std::vector<MarkerPoint>& markers)
{
    for (MarkerPoint& marker : markers) {
        for (int dir = 0; dir < 3; ++dir) {
            marker.residual[dir] = 0.0;
            marker.lambda[dir] = 0.0;
            marker.search[dir] = 0.0;
            marker.response[dir] = 0.0;
        }
    }
}

amrex::Real
markerDot(std::vector<MarkerPoint> const& markers, MarkerField a, MarkerField b)
{
    amrex::Real dot = 0.0;
    for (MarkerPoint const& marker : markers) {
        auto const& av = markerField(marker, a);
        auto const& bv = markerField(marker, b);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            dot += marker.quadrature_weight * av[dir] * bv[dir];
        }
    }
    return dot;
}

GeneralizedForce
markerGeneralizedForce(std::vector<MarkerPoint> const& markers, MarkerField field)
{
    GeneralizedForce rhs;
    for (MarkerPoint const& marker : markers) {
        auto const& values = markerField(marker, field);
        const amrex::Real weight = marker.quadrature_weight;
        const amrex::Real rx = marker.center_vector[0];
        const amrex::Real ry = marker.center_vector[1];

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            rhs.force[dir] -= weight * values[dir];
        }

#if (AMREX_SPACEDIM == 2)
        rhs.torque[2] -= weight * (rx * values[1] - ry * values[0]);
#else
        const amrex::Real rz = marker.center_vector[2];
        rhs.torque[0] -= weight * (ry * values[2] - rz * values[1]);
        rhs.torque[1] -= weight * (rz * values[0] - rx * values[2]);
        rhs.torque[2] -= weight * (rx * values[1] - ry * values[0]);
#endif
    }
    return rhs;
}

amrex::Real&
constraintField(ConstraintPoint& constraint, MarkerField field)
{
    switch (field) {
    case MarkerField::Lambda:
        return constraint.lambda;
    case MarkerField::Residual:
        return constraint.residual;
    case MarkerField::Search:
        return constraint.search;
    case MarkerField::Response:
        return constraint.response;
    }
    return constraint.residual;
}

amrex::Real const&
constraintField(ConstraintPoint const& constraint, MarkerField field)
{
    switch (field) {
    case MarkerField::Lambda:
        return constraint.lambda;
    case MarkerField::Residual:
        return constraint.residual;
    case MarkerField::Search:
        return constraint.search;
    case MarkerField::Response:
        return constraint.response;
    }
    return constraint.residual;
}

bool
hasNonzeroConstraintField(std::vector<ConstraintPoint> const& constraints, MarkerField field)
{
    for (ConstraintPoint const& constraint : constraints) {
        if (std::abs(constraintField(constraint, field)) > 0.0) {
            return true;
        }
    }
    return false;
}

void
copyConstraintField(std::vector<ConstraintPoint>& constraints, MarkerField dst, MarkerField src)
{
    for (ConstraintPoint& constraint : constraints) {
        constraintField(constraint, dst) = constraintField(constraint, src);
    }
}

void
resetConstraintUzawaFields(std::vector<ConstraintPoint>& constraints)
{
    for (ConstraintPoint& constraint : constraints) {
        constraint.residual = 0.0;
        constraint.search = 0.0;
        constraint.response = 0.0;
    }
}

amrex::Real
constraintDot(std::vector<ConstraintPoint> const& constraints, MarkerField a, MarkerField b)
{
    amrex::Real dot = 0.0;
    for (ConstraintPoint const& constraint : constraints) {
        dot += constraint.quadrature_weight
             * constraintField(constraint, a)
             * constraintField(constraint, b);
    }
    return dot;
}

GeneralizedForce
constraintGeneralizedForce(std::vector<ConstraintPoint> const& constraints, MarkerField field)
{
    return constraintGeneralizedForce(constraints, field, static_cast<std::size_t>(-1));
}

GeneralizedForce
constraintGeneralizedForce(std::vector<ConstraintPoint> const& constraints,
                           MarkerField field,
                           std::size_t body_index)
{
    GeneralizedForce rhs;
    for (ConstraintPoint const& constraint : constraints) {
        if (body_index != static_cast<std::size_t>(-1) &&
            constraint.body_index != body_index) {
            continue;
        }
        const int comp = constraint.component;
        const amrex::Real value = constraintField(constraint, field);
        const amrex::Real weight = constraint.quadrature_weight;
        const amrex::Real rx = constraint.center_vector[0];
        const amrex::Real ry = constraint.center_vector[1];

        rhs.force[comp] -= weight * value;

#if (AMREX_SPACEDIM == 2)
        if (comp == 0) {
            rhs.torque[2] += weight * ry * value;
        } else {
            rhs.torque[2] -= weight * rx * value;
        }
#else
        const amrex::Real rz = constraint.center_vector[2];
        if (comp == 0) {
            rhs.torque[1] -= weight * rz * value;
            rhs.torque[2] += weight * ry * value;
        } else if (comp == 1) {
            rhs.torque[0] += weight * rz * value;
            rhs.torque[2] -= weight * rx * value;
        } else {
            rhs.torque[0] -= weight * ry * value;
            rhs.torque[1] += weight * rx * value;
        }
#endif
    }
    return rhs;
}

void
updateLambdaAndResidual(std::vector<MarkerPoint>& markers, amrex::Real alpha)
{
    for (MarkerPoint& marker : markers) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            marker.lambda[dir] -= alpha * marker.search[dir];
            marker.residual[dir] -= alpha * marker.response[dir];
        }
    }
}

void
updateSearch(std::vector<MarkerPoint>& markers, amrex::Real beta)
{
    for (MarkerPoint& marker : markers) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            marker.search[dir] = marker.residual[dir] + beta * marker.search[dir];
        }
    }
}

void
updateLambdaAndResidual(std::vector<ConstraintPoint>& constraints, amrex::Real alpha)
{
    for (ConstraintPoint& constraint : constraints) {
        constraint.lambda -= alpha * constraint.search;
        constraint.residual -= alpha * constraint.response;
    }
}

void
updateSearch(std::vector<ConstraintPoint>& constraints, amrex::Real beta)
{
    for (ConstraintPoint& constraint : constraints) {
        constraint.search = constraint.residual + beta * constraint.search;
    }
}

void
subtractRigidResponseFromMarkers(std::vector<MarkerPoint>& markers,
                                 RigidMotion const& motion,
                                 DLMOptions const& options)
{
    for (MarkerPoint& marker : markers) {
        const auto rigid_velocity = rigidVelocityAt(marker.position, motion, options);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            marker.response[dir] -= rigid_velocity[dir];
        }
    }
}

void
subtractRigidResponseFromConstraints(std::vector<ConstraintPoint>& constraints,
                                     RigidMotion const& motion,
                                     DLMOptions const& options)
{
    for (ConstraintPoint& constraint : constraints) {
        const auto rigid_velocity = rigidVelocityAt(constraint.position, motion, options);
        constraint.response -= rigid_velocity[constraint.component];
    }
}

void
subtractRigidResponseFromConstraints(std::vector<ConstraintPoint>& constraints,
                                     std::vector<RigidMotion> const& motions,
                                     std::vector<BodyState> const& bodies)
{
    for (ConstraintPoint& constraint : constraints) {
        if (constraint.body_index >= motions.size() || constraint.body_index >= bodies.size()) {
            continue;
        }
        const auto rigid_velocity =
            rigidVelocityAt(constraint.position,
                            motions[constraint.body_index],
                            bodies[constraint.body_index]);
        constraint.response -= rigid_velocity[constraint.component];
    }
}

void
fillMarkerResiduals(std::vector<MarkerPoint>& markers,
                    RigidMotion const& motion,
                    DLMOptions const& options)
{
    for (MarkerPoint& marker : markers) {
        const auto rigid_velocity = rigidVelocityAt(marker.position, motion, options);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            marker.residual[dir] = rigid_velocity[dir] - marker.interpolated_velocity[dir];
        }
    }
}

void
fillConstraintResiduals(std::vector<ConstraintPoint>& constraints,
                        RigidMotion const& motion,
                        DLMOptions const& options)
{
    for (ConstraintPoint& constraint : constraints) {
        const auto rigid_velocity = rigidVelocityAt(constraint.position, motion, options);
        constraint.residual = rigid_velocity[constraint.component]
                            - constraint.interpolated_value;
    }
}

void
fillConstraintResiduals(std::vector<ConstraintPoint>& constraints,
                        std::vector<BodyState> const& bodies)
{
    for (ConstraintPoint& constraint : constraints) {
        if (constraint.body_index >= bodies.size()) {
            continue;
        }
        BodyState const& body = bodies[constraint.body_index];
        const auto rigid_velocity = rigidVelocityAt(constraint.position,
                                                    bodyMotion(body),
                                                    body);
        constraint.residual = rigid_velocity[constraint.component]
                            - constraint.interpolated_value;
    }
}

} // namespace ELFF::Fluid::DLM::Detail
