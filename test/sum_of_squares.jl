@testset "sum of squares in FreeGroup *-algebra" begin
      F = Groups.FreeGroup(4)
      S = [Groups.gens(F); inv.(Groups.gens(F))]

      ID = one(F)
      RADIUS = 3
      @time E_R, sizes = Groups.wlmetric_ball(S, ID, radius=2 * RADIUS)
      @test sizes == [9, 65, 457, 3201, 22409, 156865]

      b = StarAlgebras.Basis{UInt32}(E_R)

      mstr = StarAlgebras.MTable(b, size=(sizes[RADIUS], sizes[RADIUS]))

      RG = StarAlgebra(F, b, mstr)

      g, h, k, l = S[1:4]

      G = (one(RG) - RG(g))
      @test G' == one(RG) - RG(inv(g))
      @test G' * G == MA.operate_to!(zero(G), *, G', G) == 2one(RG) - RG(g) - RG(g)'
      @test star(G * G) == G' * G'

      @testset "Sums of hermitian squares" begin
            𝕀 = one(RG)

            G = (𝕀 - RG(g))
            H = (𝕀 - RG(h))
            K = (𝕀 - RG(k))
            L = (𝕀 - RG(l))
            GH = (𝕀 - RG(g * h))
            KL = (𝕀 - RG(k * l))

            X = (2𝕀 - star(RG(g)) - RG(h))
            Y = (2𝕀 - star(RG(g * h)) - RG(k))

            @test -(2𝕀 - RG(g * h) - RG(g * h)') + 2G' * G + 2H' * H == X' * X
            @test (2𝕀 - RG(g * h) - RG(g * h)') == GH' * GH
            @test -(2𝕀 - RG(g * h * k) - RG(g * h * k)') + 2GH' * GH + 2K' * K == Y' * Y
            @test -(2𝕀 - RG(g * h * k) - RG(g * h * k)') +
                  2(GH' * GH - 2G' * G - 2H' * H) +
                  4G' * G +
                  4H' * H +
                  2K' * K == Y' * Y

            @test GH' * GH - 2G' * G - 2H' * H == -X' * X
            @test -(2𝕀 - RG(g * h * k) - RG(g * h * k)') + 4G' * G + 4H' * H + 2K' * K == 2X' * X + Y' * Y

            @test GH' * GH == 2G' * G + 2H' * H - (2𝕀 - RG(g)' - RG(h))' * (2𝕀 - RG(g)' - RG(h))
            @test KL' * KL == 2K' * K + 2L' * L - (2𝕀 - RG(k)' - RG(l))' * (2𝕀 - RG(k)' - RG(l))

            @test -(2𝕀 - RG(g * h * k * l)' - RG(g * h * k * l)) + 2 * GH' * GH + 2 * KL' * KL ==
                  (2𝕀 - RG(g * h)' - RG(k * l))' * (2𝕀 - RG(g * h)' - RG(k * l))

            @test -(2𝕀 - star(RG(g * h * k * l)) - RG(g * h * k * l)) +
                  2(2G' * G + 2H' * H - (2𝕀 - RG(g)' - RG(h))' * (2𝕀 - RG(g)' - RG(h))) +
                  2(2K' * K + 2L' * L - (2𝕀 - RG(k)' - RG(l))' * (2𝕀 - RG(k)' - RG(l))) ==
                  (2𝕀 - RG(g * h)' - RG(k * l))' * (2𝕀 - RG(g * h)' - RG(k * l))

            @test -(2𝕀 - star(RG(g * h * k * l)) - RG(g * h * k * l)) +
                  2(2G' * G + 2H' * H) +
                  2(2K' * K + 2L' * L) ==
                  (2𝕀 - RG(g * h)' - RG(k * l))' * (2𝕀 - RG(g * h)' - RG(k * l)) +
                  2(2𝕀 - RG(g)' - RG(h))' * (2𝕀 - RG(g)' - RG(h)) +
                  2(2𝕀 - RG(k)' - RG(l))' * (2𝕀 - RG(k)' - RG(l))

            @test -(2𝕀 - RG(g * h * k * l)' - RG(g * h * k * l)) +
                  2(2𝕀 - RG(g * h * k)' - RG(g * h * k)) + 2L' * L ==
                  (2𝕀 - RG(g * h * k)' - RG(l))' * (2𝕀 - RG(g * h * k)' - RG(l))

            @test 2𝕀 - RG(g * h * k)' - RG(g * h * k) ==
                  2GH' * GH + 2K' * K - (2𝕀 - star(RG(g * h)) - RG(k))' * (2𝕀 - star(RG(g * h)) - RG(k))

            @test -(2𝕀 - RG(g * h * k * l)' - RG(g * h * k * l)) +
                  2(
                        2GH' * GH + 2K' * K - (2𝕀 - RG(g * h)' - RG(k))' * (2𝕀 - RG(g * h)' - RG(k))
                  ) + 2L' * L ==
                  (2𝕀 - RG(g * h * k)' - RG(l))' * (2𝕀 - RG(g * h * k)' - RG(l))

            @test -(2𝕀 - RG(g * h * k * l)' - RG(g * h * k * l)) + 2(2GH' * GH + 2K' * K) + 2L' * L ==
                  (2𝕀 - RG(g * h * k)' - RG(l))' * (2𝕀 - RG(g * h * k)' - RG(l)) +
                  2(2𝕀 - RG(g * h)' - RG(k))' * (2𝕀 - RG(g * h)' - RG(k))

            @test -(2𝕀 - RG(g * h * k * l)' - RG(g * h * k * l)) +
                  8G' * G + 8H' * H + 4K' * K + 2L' * L ==
                  (2𝕀 - RG(g * h * k)' - RG(l))' * (2𝕀 - RG(g * h * k)' - RG(l)) +
                  2(2𝕀 - RG(g * h)' - RG(k))' * (2𝕀 - RG(g * h)' - RG(k)) +
                  4(2𝕀 - RG(g)' - RG(h))' * (2𝕀 - RG(g)' - RG(h))

            @test -(2𝕀 - RG(g * h * k * l)' - RG(g * h * k * l)) + 2GH' * GH + 2KL' * KL ==
                  (2𝕀 - RG(g * h)' - RG(k * l))' * (2𝕀 - RG(g * h)' - RG(k * l))

            @test -(2𝕀 - RG(g * h * k * l)' - RG(g * h * k * l)) +
                  2(2G' * G + 2H' * H) +
                  2(2K' * K + 2L' * L) ==
                  (2𝕀 - RG(g * h)' - RG(k * l))' * (2𝕀 - RG(g * h)' - RG(k * l)) +
                  2(2𝕀 - RG(k)' - RG(l))' * (2𝕀 - RG(k)' - RG(l)) +
                  2(2𝕀 - RG(g)' - RG(h))' * (2𝕀 - RG(g)' - RG(h))
      end
end
