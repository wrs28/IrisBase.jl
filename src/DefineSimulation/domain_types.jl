export AbstractDomain,
GenericDomain

export isWaveguide,
isBulkWaveguide,
isBackground,
isDefect,
isPC,
isPlanar,
isHalfSpace

abstract type AbstractForegroundBackgound end
struct Foreground<:AbstractForegroundBackgound end
struct Background<:AbstractForegroundBackgound end

abstract type AbstractWaveguideNot end
struct Waveguide<:AbstractWaveguideNot end
struct NotWaveguide<:AbstractWaveguideNot end

abstract type AbstractDomainType end
struct PlanarDomainType<:AbstractDomainType end
struct HalfspaceDomainType<:AbstractDomainType end
struct PCDomainType<:AbstractDomainType end
struct DefectDomainType<:AbstractDomainType end
struct GenericDomainType<:AbstractDomainType end

abstract type AbstractBulkAsymptote end
struct Bulk<:AbstractBulkAsymptote end
struct Asymptote<:AbstractBulkAsymptote end


abstract type AbstractDomain{TFB,TWN,TDT,TBA} end
struct PlanarWaveguide          <:AbstractDomain{Foreground,Waveguide   ,PlanarDomainType   ,Asymptote  } end
Base.show(io::IO,::PlanarWaveguide) = print("Planar Waveguide")
struct PlanarWaveguideBackground<:AbstractDomain{Background,Waveguide   ,PlanarDomainType   ,Asymptote  } end
Base.show(io::IO,::PlanarWaveguideBackground) = print("Planar Waveguide Background")
struct Halfspace                <:AbstractDomain{Foreground,Waveguide   ,HalfspaceDomainType,Asymptote  } end
Base.show(io::IO,::Halfspace) = print("Halfspace")
struct PCWaveguide              <:AbstractDomain{Foreground,Waveguide   ,PCDomainType       ,Asymptote  } end
Base.show(io::IO,::PCWaveguide) = print("PC Waveguides")
struct PCWaveguideBackground    <:AbstractDomain{Background,Waveguide   ,PCDomainType       ,Asymptote  } end
Base.show(io::IO,::PCWaveguideBackground) = print("PC Waveguide Background")
struct BulkPlanarWaveguideX     <:AbstractDomain{Foreground,Waveguide   ,PlanarDomainType   ,Bulk       } end
Base.show(io::IO,::BulkPlanarWaveguideX) = print("Bulk Planar Waveguide X")
struct BulkPlanarWaveguideY     <:AbstractDomain{Foreground,Waveguide   ,PlanarDomainType   ,Bulk       } end
Base.show(io::IO,::BulkPlanarWaveguideY) = print("Bulk Planar Waveguide Y")
struct BulkPCWaveguideX         <:AbstractDomain{Foreground,Waveguide   ,PCDomainType       ,Bulk       } end
Base.show(io::IO,::BulkPCWaveguideX) = print("Bulk PC Waveguide X")
struct BulkPCWaveguideY         <:AbstractDomain{Foreground,Waveguide   ,PCDomainType       ,Bulk       } end
Base.show(io::IO,::BulkPCWaveguideY) = print("Bulk PC Waveguide Y")
struct BackgroundDomain         <:AbstractDomain{Background,NotWaveguide,GenericDomainType  ,Bulk       } end
Base.show(io::IO,::BackgroundDomain) = print("Background")
struct DefectDomain             <:AbstractDomain{Foreground,NotWaveguide,DefectDomainType   ,Bulk       } end
Base.show(io::IO,::DefectDomain) = print("Defect")
struct SiteDefect               <:AbstractDomain{Foreground,NotWaveguide,DefectDomainType   ,Bulk       } end
Base.show(io::IO,::SiteDefect) = print("Site Defect")
struct LineDefect               <:AbstractDomain{Foreground,NotWaveguide,DefectDomainType   ,Bulk       } end
Base.show(io::IO,::LineDefect) = print("Line Defect")
struct PCDomain                 <:AbstractDomain{Foreground,NotWaveguide,PCDomainType       ,Bulk       } end
Base.show(io::IO,::PCDomain) = print("PC")
struct GenericDomain            <:AbstractDomain{Background,NotWaveguide,GenericDomainType  ,Bulk       } end
Base.show(io::IO,::GenericDomain) = print("Generic")


isWaveguide(domain::T) where T<:AbstractDomain{_1,Waveguide} where _1= true
isWaveguide(domain) = false

isBulkWaveguide(domain::T) where T<:AbstractDomain{_1,Waveguide,_2,Bulk} where {_1,_2} = true
isBulkWaveguide(domain) = false

isBackground(domain::T) where T<:AbstractDomain{Background} = true
isBackground(domain) = false

isDefect(domain::T) where T<:AbstractDomain{_1,_2,DefectDomainType} where {_1,_2} = true
isDefect(domain) = false

isPC(domain::T) where T<:AbstractDomain{_1,_2,PCDomainType} where {_1,_2} = true
isPC(domain) = false

isPlanar(domain::T) where T<:AbstractDomain{_1,_2,PlanarDomainType} where {_1,_2} = true
isPlanar(domain) = false

isHalfSpace(domain::T) where T<:AbstractDomain{_1,_2,HalfspaceDomainType} where {_1,_2} = true
isHalfSpace(domain) = false
