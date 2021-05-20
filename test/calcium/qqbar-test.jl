@testset "qqbar.constructors" begin
   R = CalciumQQBar

   @test elem_type(R) == qqbar
   @test elem_type(CalciumQQBarField) == qqbar
   @test parent_type(qqbar) == CalciumQQBarField

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
end


@testset "qqbar.manipulation" begin
   R = CalciumQQBar

   @test iszero(R(0))
   @test isone(R(1))
   @test isrational(R(1))
   @test isreal(R(1))
   @test degree(R(1)) == 1

   u = sqrt(R(2))
   i = sqrt(R(-1))

   @test degree(u) == 2
   @test !isrational(u)
   @test isreal(u)
   @test !isrational(i)
   @test !isreal(i)
   @test is_algebraic_integer(u)

   ZZx, x = PolynomialRing(FlintZZ, "x")
   QQy, y = PolynomialRing(FlintQQ, "x")

   @test minpoly(ZZx, u) == x^2 - 2
   @test minpoly(QQy, u) == y^2 - 2

   @test root(qqbar(-1), 3) == root_of_unity(R, 6)
   @test root_of_unity(R, 4) == i
   @test root_of_unity(R, 4, 3) == -i

   @test_throws DivideError (R(1) // R(0))
   @test_throws DomainError (R(0) ^ R(-1))
   @test_throws DomainError (root(R(1), 0))
   @test_throws DomainError (u ^ u)
   @test_throws DomainError (root_of_unity(R, 0))
   @test_throws DomainError (root_of_unity(R, 0, 1))

   v = roots(x^5-x-1, CalciumQQBar)
   @test v[1]^5 - v[1] - 1 == 0

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


