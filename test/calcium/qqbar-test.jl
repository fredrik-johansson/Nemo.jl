@testset "qqbar.constructors" begin
   R = CalciumQQBar

   @test elem_type(R) == qqbar
   @test elem_type(CalciumQQBarField) == qqbar
   @test parent_type(qqbar) == CalciumQQBarField
   @test isdomain_type(qqbar) == true
   @test base_ring(CalciumQQBar) == CalciumQQBar
   @test base_ring(qqbar(3)) == CalciumQQBar

   @test isa(R, CalciumQQBarField)

   @test isa(R(), qqbar)
   @test isa(R(2), qqbar)
   @test isa(R(fmpz(2)), qqbar)
   @test isa(R(fmpq(2)), qqbar)

   @test isa(qqbar(), qqbar)
   @test isa(qqbar(2), qqbar)
   @test isa(qqbar(fmpz(2)), qqbar)
   @test isa(qqbar(fmpq(2)), qqbar)

end

@testset "qqbar.printing" begin
   a = CalciumQQBar(1)

   @test string(a) == "1.00000 (deg 1)"
   @test string(parent(a)) == "Field of Algebraic Numbers in minimal polynomial representation"
   @test needs_parentheses(a) == true
end


@testset "qqbar.manipulation" begin
   R = CalciumQQBar

   @test zero(R) == 0
   @test one(R) == 1
   @test isa(zero(R), qqbar)
   @test isa(one(R), qqbar)

   @test iszero(R(0))
   @test isone(R(1))
   @test isrational(R(1))
   @test isreal(R(1))
   @test degree(R(1)) == 1

   u = sqrt(R(2))
   i = sqrt(R(-1))

   @test canonical_unit(u) == u

   @test degree(u) == 2
   @test !isrational(u)
   @test isreal(u)
   @test !isrational(i)
   @test !isreal(i)
   @test is_algebraic_integer(u)

   @test abs(-u) == u
   @test abs2(u) == 2
   @test u != i
   @test sign(2*i) == i
   @test conj(i) == -i
   @test real(3+4*i) == 3
   @test imag(3+4*i) == 4
   @test csgn(i) == 1
   @test sign_real(-3+4*i) == -1
   @test sign_imag(-3+4*i) == 1
   @test floor(u) == 1
   @test ceil(u) == 2

   @test (u >> 3) == u // 8
   @test (u << 3) == 8 * u

   @test u < 2
   @test u > 1
   @test_throws DomainError (i > 1)

   ZZx, x = PolynomialRing(FlintZZ, "x")
   QQy, y = PolynomialRing(FlintQQ, "x")

   @test minpoly(ZZx, u) == x^2 - 2
   @test minpoly(QQy, u) == y^2 - 2

   @test root(qqbar(-1), 3) == root_of_unity(R, 6)
   @test root_of_unity(R, 4) == i
   @test root_of_unity(R, 4, 3) == -i

   @test sinpi(qqbar(1)//6) == qqbar(1)//2
   @test cospi(qqbar(1)//3) == qqbar(1)//2
   @test tanpi(qqbar(1)//3) == sqrt(qqbar(3))
   @test_throws DomainError tanpi(qqbar(1)//2)

   @test atanpi(sqrt(qqbar(3))) == qqbar(1)//3
   @test asinpi(sqrt(qqbar(2))//2) == qqbar(1)//4
   @test acospi(sqrt(qqbar(3))//2) == qqbar(1)//6
   @test log_pi_i(i) == qqbar(1)//2

   @test_throws DomainError atanpi(qqbar(2))
   @test_throws DomainError asinpi(qqbar(2))
   @test_throws DomainError acospi(qqbar(2))
   @test_throws DomainError log_pi_i(qqbar(2))

   @test_throws DivideError (R(1) // R(0))
   @test_throws DomainError (R(0) ^ R(-1))
   @test_throws DomainError (root(R(1), 0))
   @test_throws DomainError (u ^ u)
   @test_throws DomainError (root_of_unity(R, 0))
   @test_throws DomainError (root_of_unity(R, 0, 1))

   @test is_root_of_unity(i)
   @test !is_root_of_unity(qqbar(2))
   @test root_of_unity_as_args(-i) == (4, 3)
   @test_throws DomainError root_of_unity_as_args(qqbar(2))

   v = roots(x^5-x-1, CalciumQQBar)
   @test v[1]^5 - v[1] - 1 == 0

   v = roots(y^2+1, CalciumQQBar)
   @test v == [i, -i]

   @test roots(ZZx(0), CalciumQQBar) == []
   @test roots(ZZx(1), CalciumQQBar) == []
   @test roots(QQy(0), CalciumQQBar) == []
   @test roots(QQy(1), CalciumQQBar) == []

   @test eigenvalues(zero(MatrixSpace(ZZ, 0, 0)), CalciumQQBar) == []
   @test eigenvalues(zero(MatrixSpace(QQ, 0, 0)), CalciumQQBar) == []
   @test eigenvalues(ZZ[1 1; 1 -1], CalciumQQBar) == [u, -u]
   @test eigenvalues(QQ[1 1; 1 -1], CalciumQQBar) == [u, -u]

   @test conjugates(qqbar(3)) == [qqbar(3)]
   @test conjugates(u) == [u, -u]

   @test fmpz(qqbar(3)) == 3
   @test fmpq(qqbar(3) // 2) == fmpq(3,2)

   RR = ArbField(128)
   CC = AcbField(128)

   @test RR(qqbar(3)) == 3
   @test CC(qqbar(3)) == 3
   @test_throws DomainError (RR(i))

   v = sqrt(RR(2)) + sqrt(RR(3))
   @test guess(CalciumQQBar, v, 4) == sqrt(qqbar(2)) + sqrt(qqbar(3))
   @test guess(CalciumQQBar, v, 4, 10) == sqrt(qqbar(2)) + sqrt(qqbar(3))
   @test_throws ErrorException guess(CalciumQQBar, v, 2)

   @test guess(CalciumQQBar, CC(2+i), 2, 10) == 2+i

end

@testset "qqbar.adhoc-operations" begin
   R = CalciumQQBar

   @test qqbar(2) + qqbar(3) == 5
   @test qqbar(2) + 3 == 5
   @test qqbar(2) + fmpz(3) == 5
   @test qqbar(2) + fmpq(3) == 5
   @test 3 + qqbar(2) == 5
   @test fmpz(3) + qqbar(2) == 5
   @test fmpq(3) + qqbar(2) == 5

   @test qqbar(2) - qqbar(3) == -1
   @test qqbar(2) - 3 == -1
   @test qqbar(2) - fmpz(3) == -1
   @test qqbar(2) - fmpq(3) == -1
   @test 3 - qqbar(2) == 1
   @test fmpz(3) - qqbar(2) == 1
   @test fmpq(3) - qqbar(2) == 1

   @test qqbar(2) * qqbar(3) == 6
   @test qqbar(2) * 3 == 6
   @test qqbar(2) * fmpz(3) == 6
   @test qqbar(2) * fmpq(3) == 6
   @test 3 * qqbar(2) == 6
   @test fmpz(3) * qqbar(2) == 6
   @test fmpq(3) * qqbar(2) == 6

   @test qqbar(6) // qqbar(2) == 3
   @test qqbar(6) // 2 == 3
   @test qqbar(6) // fmpz(2) == 3
   @test qqbar(6) // fmpq(2) == 3
   @test 6 // qqbar(2) == 3
   @test fmpz(6) // qqbar(2) == 3
   @test fmpq(6) // qqbar(2) == 3

   @test qqbar(2) ^ qqbar(3) == 8
   @test qqbar(2) ^ 3 == 8
   @test qqbar(2) ^ fmpz(3) == 8
   @test qqbar(2) ^ fmpq(3) == 8
   @test 2 ^ qqbar(3) == 8
   @test fmpz(2) ^ qqbar(3) == 8
   @test fmpq(2) ^ qqbar(3) == 8

   @test qqbar(2) < qqbar(3)
   @test qqbar(2) < 3
   @test qqbar(2) < fmpz(3)
   @test qqbar(2) < fmpq(3)
   @test 2 < qqbar(3)
   @test fmpz(2) < qqbar(3)
   @test fmpq(2) < qqbar(3)

end

@testset "qqbar.rand" begin
   R = CalciumQQBar

   for i=1:10
      x = rand(CalciumQQBar, degree=5, bits=5)
      @test degree(x) <= 5
   end

   for i=1:10
      x = rand(CalciumQQBar, degree=5, bits=5, randtype=:real)
      @test isreal(x)
   end

   for i=1:10
      x = rand(CalciumQQBar, degree=5, bits=5, randtype=:nonreal)
      # todo: need to upgrade Calcium
      # @test !isreal(x)
   end

end


