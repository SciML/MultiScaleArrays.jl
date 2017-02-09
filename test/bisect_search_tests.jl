a = collect(3:3:30)
@test MultiScaleArrays.bisect_search(a,20) == 7
@test MultiScaleArrays.bisect_search(a,13) == 5
@test MultiScaleArrays.bisect_search(a,12) == 4
for i = 1:30
  @test MultiScaleArrays.bisect_search(a,i) == ((i-1)÷3)+1 #+1 for 1-based indexing
end
