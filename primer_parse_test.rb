a = Read.new("test_name", "+", true, false, "This is ANOTHER example read", [])
puts a.inspect

b= Primer.new(1, 100, 2, "+", "todo")
puts b.inspect

a.primers.push(b)

puts "After adding primer:"
puts a.inspect
