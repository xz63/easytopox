function area = triangle(a, b, c)
s = (a+b+c)/2;
area = sqrt(s*(s-a)*(s-b)*(s-c));
