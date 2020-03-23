test_fun = function(a, b, c){

  a + b + c

}


test_list = list("a" = 1, "b" = 3)

do.call(test_fun, c(test_list, "c" = 100))
