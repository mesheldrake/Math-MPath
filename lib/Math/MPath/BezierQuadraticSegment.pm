####################################################################################
###      Math::MPath::BezierQuadraticSegment   ######################################
package Math::MPath::BezierQuadraticSegment;
{
# until we do this, maybe we get away with using cubic Beziers with both
# tangent handle points on the same spot.
push @Math::MPath::BezierQuadraticSegment::ISA, 'Math::MPath::BezierCubicSegment';
}
1;