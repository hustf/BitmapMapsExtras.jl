
z = z_cylinder(0.1)[38:end, 106:end]

@testset "NegateY" begin
    # 6 high, 10 wide matrix
    R = CartesianIndices((6, 10))
    flipy = NegateY(R)
    @test flipy(0) == 7
    @test flipy(1) == 6
    @test flipy(2) == 5
    @test flipy(6) == 1
    @test flipy(6.0) == 1
end

@testset "NegateY reversible" begin
    R = CartesianIndices((6, 10))
    flipy = NegateY(R)
    y = 1.0
    i = flipy(y)
    @test flipy(i) == y
    y = 6.0
    i = flipy(y)
    @test flipy(i) == y
end

@testset "Domain" begin
    # 8 high, 10 wide matrix
    R = CartesianIndices((8, 10))
    d = Domain(R)
    @test d == Domain(1.0, 1.0, 10.0, 8.0)
    #
    Ω = CartesianIndices((-2:2, -2:2))
    d = Domain(R, Ω)
    @test d == Domain(3.0, 3.0, 8.0, 6.0)
    @test !d(0.5, 0.3)
    @test !d(0.5, 3)
    @test d(5, 3)
    @test d(5.5, 3.5)
    @test !d(1000, 32)
    @test !d(30, 3200)
    @test !d(3045, 3200)
    #
    @test signed_distance_within_domain(d, 3.0, 4.0) == 0.0
    @test signed_distance_within_domain(d, 8.0, 4.0) == 0.0
    @test signed_distance_within_domain(d, 4.0, 3.0) == 0.0
    @test signed_distance_within_domain(d, 4.0, 6.0) == 0.0
    #
    @test signed_distance_within_domain(d, 3.5, 4.0) == 0.5
    @test signed_distance_within_domain(d, 7.5, 4.0) == 0.5
    @test signed_distance_within_domain(d, 4.0, 3.5) == 0.5
    @test signed_distance_within_domain(d, 4.0, 5.5) == 0.5
    #
    @test signed_distance_within_domain(d, 2.5, 4.0) == -0.5
    @test signed_distance_within_domain(d, 8.5, 4.0) == -0.5
    @test signed_distance_within_domain(d, 4.0, 2.5) == -0.5
    @test signed_distance_within_domain(d, 4.0, 6.5) == -0.5
end




@testset "𝐧ₚ! Surface normal vector projection" begin 
    M = z[11:15, 11:15]
    v = MVector{2, Float64}([0, 0])
    𝐧ₚ!(v, M)
    @test v[1] < 0 # Normal vector points left
    @test v[2] > 0 # Normal vector points up ("y is up")
    @test abs(v[1]) < abs(v[2]) # Small horizontal component, large vertical
end

@testset "DirectionOnGrid" begin
    dog = DirectionOnGrid(𝐧ₚ!, z)
    z[5,5]
    v = MVector{2, Float64}([0, 0])
    pt = CartesianIndex(5,5)
    direction_on_grid!(v, dog, pt)
    @test v[1] < 0 # Normal vector points left
    @test v[2] > 0 # Normal vector points up ("y is up")
    @test v ≈ [-0.09903836842989576, 0.9865965311601863]
    pt = CartesianIndex(5,6)
    direction_on_grid!(v, dog, pt)
    @test v ≈ [-0.09901472874636574, 0.9863907566538036]
end

@testset "DirectionInDomain internals" begin 
    did = DirectionInDomain(𝐧ₚ!, z)
    did.li.coefs[2,2] .= [1.0, 10.0]
    did.li(1.5, 1.5) == [0.25, 2.5]
    did.li(1.0, 1.0) == [0., 0.]
    did.li(2.0, 2.0) == [1.0, 10.]
    typeof(did.li(2.0, 2.0)) <: MVector
end
@testset "DirectionInDomain" begin
    did = DirectionInDomain(𝐧ₚ!, z)
    @test DirectionInDomain(𝐧ₚ!, z)(5.0, 5.0) ≈ [-0.09903836842989576, 0.9865965311601863]
    v = MVector{2, Float64}([0, 0])
    direction_in_domain!(v, did, 5.0, 5.0)
    @test v ≈ [-0.09903836842989576, 0.9865965311601863]
    direction_in_domain!(v, did, 6.0, 5.0)
    @test v ≈ [-0.09901472874636574, 0.9863907566538036]
end

@testset "DirectionAtXY" begin
    daxy = DirectionAtXY(𝐧ₚ!, z)
    v = daxy(3.0, 3.0)
    @test v ≈ [0.09110987888453688, -0.9080563695155759]
    @test_throws AssertionError daxy(2.99, 3.0)
    v = daxy(892.0, 960.0)
    @test v ≈ [-0.0816455901667088, 0.8137317920512829]
    @test_throws AssertionError daxy(892.1, 960.0)
    @test daxy.d == Domain(3.0, 3.0, 892.0, 960.0)
    # This point up and a little to the left
    𝐧ₚx, 𝐧ₚy = daxy(962 / 2, 894 - 100)
    @test abs(𝐧ₚx) < abs(𝐧ₚy)
    @test 𝐧ₚx < 0 # Normal vector points left
    @test 𝐧ₚy > 0 # Normal vector points up ("y is up")
end

# daxy = DirectionAtXY(𝐧ₚ!, z)
# 269.040 ns (9 allocations: 336 bytes)
# @btime daxy(892.0, 960.0)

# Visualize points, no traces
function test_plot(fdir!, z, Δ, lenmult)
    # Background: z-values in color
    img = bluesc(z; mi = -float(r) ) .* 1.0
    # Add simple contour lines, too
    Δc = 20
    wc = Δc / 10
    map!(img, z, img) do zz, pix
        mod(zz, Δc) < wc ? RGBA{N0f8}(0.1, 0.1, 0.1, 1.0) : pix 
    end
    test_plot!(img, fdir!,  z, Δ, lenmult)
end
function test_plot!(img, fdir!, z, Δ, lenmult)
    # Black-white buffer
    bbuf = Array{GrayA{N0f8}}(falses( size(img)...))
    test_plot!(bbuf, fdir!, z, Δ, lenmult)
    # Function converting GrayA{N0f8} to proper color
    f = x -> RGBA{N0f8}(COLOR_CURVGLYPH.r, COLOR_CURVGLYPH.g, COLOR_CURVGLYPH.b, 
                            x.val)
    # Composite bbuf over img
    # Overlay bbuf on img in the proper color
    map!(BlendLighten, img, img, f.(bbuf))
    img
end

function test_plot!(bbuf::Array{GrayA{N0f8}}, fdir!, z, Δ, lenmult)
    daxy = DirectionAtXY(fdir!, z)
    d = daxy.d
    ny = daxy.ny
    for x in d.minx:Δ:d.maxx, y in d.miny:Δ:d.maxy
        𝐧ₚx, 𝐧ₚy = daxy(x, y)
        Δi = Int(round(lenmult * -𝐧ₚy))
        Δj = Int(round(lenmult * 𝐧ₚx))
        i = Int(round(ny(y)))
        j = Int(round(x))
        pt = CartesianIndex(i, j)
        vector!(bbuf, pt, Δi, Δj)
    end
    bbuf
end

test_plot(𝐧ₚ!, z, 50, 50 )
test_plot(𝐧ₚ!, z_sphere(), 50, 50 )
test_plot(𝐧ₚ!, z_ellipsoid(), 50, 50 )
test_plot(𝐧ₚ!, z_paraboloid(; a = r, b= 0.5r), 50, 50 )
test_plot(𝐧ₚ!, z_paraboloid(; a = -r, b= -0.5r), 50, 50 )
test_plot(𝐧ₚ!, z_paraboloid(; a = -2r, b= 2r), 50, 150 )
test_plot(𝐧ₚ!, z_paraboloid(; a = 2r, b= -2r), 50, 150 )


