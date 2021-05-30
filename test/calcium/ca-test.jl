@testset "ca.constructors" begin
   C = CalciumField()

   @test elem_type(C) == ca
   @test elem_type(CalciumField) == ca
   @test parent_type(ca) == CalciumField
   @test isdomain_type(ca) == true
   @test base_ring(C) == Union{}      # ?
   @test base_ring(C(3)) == Union{}      # ?

   @test isa(C, CalciumField)

   @test isa(C(), ca)
   @test isa(C(2), ca)
   @test isa(C(2+3im), ca)
   @test isa(C(fmpz(2)), ca)
   @test isa(C(fmpq(2)), ca)
   @test isa(C(qqbar(2)), ca)
   @test isa(C(C(2)), ca)

   C2 = CalciumField()

   a = C(3)
   a2 = C2(3)

   @test parent(a) == C
   @test parent(a2) == C2
   @test parent(parent(a)(a2)) == C
   @test parent(parent(a2)(a)) == C2

end

@testset "ca.printing" begin
   C = CalciumField()
   Cext = CalciumField(extended=true)

   @test string(C) == "Exact Complex Field"
   @test string(Cext) == "Exact Complex Field (Extended)"

   @test string(C(1)) == "1"
   @test string(Cext(1)) == "1"

   @test string(C(pi)) == "3.14159 {a where a = 3.14159 [Pi]}"

end

@testset "ca.manipulation" begin
   C = CalciumField()
   Cext = CalciumField(extended=true)

   @test zero(C) == 0
   @test one(C) == 1
   @test isa(zero(C), ca)
   @test isa(one(C), ca)

   @test iszero(C(0))
   @test isone(C(1))
   @test isinteger(C(1))
   @test isrational(C(1))
   @test isreal(C(1))
   @test isnumber(C(1))

   u = sqrt(C(2))
   i = sqrt(C(-1))

   @test i == C(0+1im)
   @test 3+4*i == C(3+4im)

   @test u == Cext(u)
   @test u == C(Cext(u))
   u_i = u + i
   @test u_i == Cext(u_i)
   @test u_i == C(Cext(u_i))

   @test canonical_unit(u) == u
   @test isa(hash(u), UInt)

   @test !isinteger(u)
   @test !isrational(u)
   @test isreal(u)
   @test !isrational(i)
   @test !isreal(i)
   @test isimaginary(i)
   @test !isimaginary(u)

   @test inv(u) == u // 2

   @test abs(-u) == u
   @test u != i
   @test sign(2*i) == i
   @test conj(i) == -i
   @test real(3+4*i) == 3
   @test imag(3+4*i) == 4
   @test csgn(i) == 1
   #@test sign_real(-3+4*i) == -1
   #@test sign_imag(-3+4*i) == 1
   @test floor(u) == 1
   @test ceil(u) == 2

   @test_throws DomainError infinity(C)
   @test_throws DomainError unsigned_infinity(C)
   @test_throws DomainError undefined(C)
   @test_throws DomainError unknown(C)

   inf = infinity(Cext)
   uinf = unsigned_infinity(Cext)
   und = undefined(Cext)
   unk = unknown(Cext)

   @test_throws DomainError C(inf)

   @test_throws DomainError C(1) // 0
   @test Cext(1) // 0 == uinf

   @test_throws DomainError log(C(0))
   @test log(Cext(0)) == -inf

   @test_throws DomainError C(0) // 0
   @test Cext(0) // 0 == undefined(Cext)

   @test -2*Cext(i)*inf == infinity(Cext(-i))

   @test is_signed_inf(inf)
   @test !is_signed_inf(uinf)
   @test isinf(inf)
   @test isinf(uinf)
   @test !isuinf(inf)
   @test isuinf(uinf)

   @test isundefined(und)
   @test !isunknown(und)

   @test isunknown(unk)
   @test_throws ErrorException isreal(unk)
   @test_throws ErrorException isnumber(unk)
   @test_throws ErrorException isundefined(unk)

   @test !isinf(C(1))
   @test !isuinf(C(1))
   @test !is_signed_inf(C(1))
   @test !isundefined(C(1))
   @test !isunknown(C(1))

   @test und == und
   @test_throws ErrorException (unk == unk)

end

@testset "ca.adhoc-operations" begin
   C = CalciumField()

   @test C(2) + C(3) == 5
   @test C(2) + 3 == 5
   @test C(2) + fmpz(3) == 5
   @test C(2) + fmpq(3) == 5
   @test C(2) + qqbar(3) == 5
   @test 3 + C(2) == 5
   @test fmpz(3) + C(2) == 5
   @test fmpq(3) + C(2) == 5
   @test qqbar(3) + C(2) == 5

   @test C(2) - C(3) == -1
   @test C(2) - 3 == -1
   @test C(2) - fmpz(3) == -1
   @test C(2) - fmpq(3) == -1
   @test C(2) - qqbar(3) == -1
   @test 3 - C(2) == 1
   @test fmpz(3) - C(2) == 1
   @test fmpq(3) - C(2) == 1
   @test qqbar(3) - C(2) == 1

   @test C(2) * C(3) == 6
   @test C(2) * 3 == 6
   @test C(2) * fmpz(3) == 6
   @test C(2) * fmpq(3) == 6
   @test C(2) * qqbar(3) == 6
   @test 3 * C(2) == 6
   @test fmpz(3) * C(2) == 6
   @test fmpq(3) * C(2) == 6
   @test qqbar(3) * C(2) == 6

   @test C(6) // C(2) == 3
   @test C(6) // 2 == 3
   @test C(6) // fmpz(2) == 3
   @test C(6) // fmpq(2) == 3
   @test C(6) // qqbar(2) == 3
   @test 6 // C(2) == 3
   @test fmpz(6) // C(2) == 3
   @test fmpq(6) // C(2) == 3
   @test qqbar(6) // C(2) == 3

   @test C(2) ^ C(3) == 8
   @test C(2) ^ 3 == 8
   @test C(2) ^ fmpz(3) == 8
   @test C(2) ^ fmpq(3) == 8
   @test 2 ^ C(3) == 8
   @test fmpz(2) ^ C(3) == 8
   @test fmpq(2) ^ C(3) == 8
   @test qqbar(2) ^ C(3) == 8

   @test C(2) < C(3)
   @test C(2) < 3
   @test C(2) < fmpz(3)
   @test C(2) < fmpq(3)
   @test C(2) < qqbar(3)
   @test 2 < C(3)
   @test fmpz(2) < C(3)
   @test fmpq(2) < C(3)
   @test qqbar(2) < C(3)

end

@testset "ca.inplace" begin
   C = CalciumField()
   C2 = CalciumField()

   x = C(7)
   zero!(x)
   @test x == 0

   x = C(7)
   y = mul!(x, C(3), C(5))
   @test x == 15
   @test x === y

   @test_throws ErrorException mul!(x, C2(3), C(5))
   @test_throws ErrorException mul!(x, C(3), C2(5))

   x = C(7)
   y = addeq!(x, C(3))
   @test x == 10
   @test x === y

   @test_throws ErrorException addeq!(x, C2(3))

   x = C(7)
   y = add!(x, C(3), C(5))
   @test x == 8
   @test x === y

   @test_throws ErrorException add!(x, C2(3), C(5))
   @test_throws ErrorException add!(x, C(3), C2(5))

end

@testset "ca.functions" begin
   C = CalciumField()

   u = sqrt(C(2))
   i = sqrt(C(-1))

   @test real(3+4*i) == 3
   @test imag(3+4*i) == 4
   @test angle(2+2*i) == C(pi) // 4
   @test csgn(-i) == -1
   @test sign(2*i) == i
   @test abs(1+i) == u
   @test conj(1+i) == 1-i
   @test conj(1+C(pi)*i, form=:deep) == 1-C(pi)*i
   @test conj(1+C(pi)*i, form=:shallow) == 1-C(pi)*i
   @test_throws ErrorException conj(1+C(pi)*i, form=:gollum)

   @test floor(u) == 1
   @test ceil(u) == 2
   @test sqrt(i) == (1+i)*u//2
   @test exp(C(pi) * i) == -1
   @test log(exp(u)) == u

   @test pow(1 + C(pi), 25) == pow(1 + C(pi), 25, form=:arithmetic)
   @test_throws ErrorException pow(1 + C(pi), 25, form=:gollum)

   @test sin(C(1)) == sin(C(1), form=:exponential)
   @test sin(C(1)) == sin(C(1), form=:tangent)
   @test sin(C(1)) == sin(C(1), form=:direct)
   @test_throws ErrorException sin(C(1), form=:gollum)

   @test cos(C(1)) == cos(C(1), form=:exponential)
   @test cos(C(1)) == cos(C(1), form=:tangent)
   @test cos(C(1)) == cos(C(1), form=:direct)
   @test_throws ErrorException cos(C(1), form=:gollum)

   @test cos(u)^2 + sin(u)^2 == 1

   @test tan(C(1)) == tan(C(1), form=:exponential)
   @test tan(C(1)) == tan(C(1), form=:sine_cosine)
   @test tan(C(1)) == tan(C(1), form=:direct)
   @test_throws ErrorException tan(C(1), form=:gollum)

   @test atan(C(1)) == C(pi)//4
   @test atan(C(2)) == atan(C(2), form=:logarithm)
   @test atan(C(2)) == atan(C(2), form=:arctangent)
   @test atan(C(2)) == atan(C(2), form=:direct)
   @test_throws ErrorException atan(C(2), form=:gollum)

   @test asin(C(1)) == C(pi)//2
   @test asin(C(2)) == asin(C(2), form=:logarithm)
   @test asin(C(2)) == asin(C(2), form=:direct)
   @test_throws ErrorException asin(C(2), form=:gollum)

   @test acos(C(-1)) == C(pi)
   @test acos(C(2)) == acos(C(2), form=:logarithm)
   @test acos(C(2)) == acos(C(2), form=:direct)
   @test_throws ErrorException acos(C(2), form=:gollum)

   @test gamma(C(5)) == 24
   @test erf(C(1)) == 1 - erfc(C(1))
   @test erfi(C(1)) == -i*erf(i)

end

