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

