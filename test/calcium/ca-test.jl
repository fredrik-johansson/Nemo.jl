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


