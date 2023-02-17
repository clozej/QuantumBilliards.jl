
abstract type CoordinateSystem end
abstract type AbsBilliard end
#abstract type PolygonalBilliard <: AbsBilliard end
abstract type AbsCurve end
abstract type AbsVirtualCurve <: AbsCurve end
abstract type AbsRealCurve <: AbsCurve end
abstract type AbsBasis end
abstract type AbsSolver end
abstract type AbsPoints end
abstract type SweepSolver <: AbsSolver end
abstract type AcceleratedSolver <: AbsSolver end
#abstract type AbsSampler end
abstract type AbsObservable end
abstract type AbsState end
abstract type StationaryState <: AbsState end
#grid iterators
abstract type AbsGrid end


