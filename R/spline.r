.intt.params = function (cy1, cy2, ct1, ct2)
{	A = cy1
	B = ct1
	C = -3 * cy1 + 3 * cy2 - 2 * ct1 - ct2
	D = 2 * cy1 - 2 * cy2 + ct1 + ct2
	c (A, B, C, D)
}

.intt = function (cx1, cx2, cy1, cy2, ct1, ct2, x)
{   if (length (x) != 1)
		stop (".intt needs scalars")
	if (!missing (cx1) && !missing (cx2) )
    {   dx = cx2 - cx1
        ct1 = dx * ct1
        ct2 = dx * ct2
        x = (x - cx1) / dx
    }
    params = .intt.params (cy1, cy2, ct1, ct2)
	sum (params * x ^ (0:3) )
}

.intt.argmax = function (cx1, cx2, cy1, cy2, ct1, ct2)
{   has.cx = (!missing (cx1) && !missing (cx2) )
	if (has.cx)
	{   dx = cx2 - cx1
        ct1 = dx * ct1
        ct2 = dx * ct2
    }

	A = 6 * cy1 - 6 * cy2 + 3 * ct1 + 3 * ct2
	B = -6 * cy1 + 6 * cy2 - 4 * ct1 - 2 * ct2
	C = ct1

	x1 = (-B + sqrt (B ^ 2 - 4 * A * C) ) / (2 * A)
	x2 = (-B - sqrt (B ^ 2 - 4 * A * C) ) / (2 * A)
	if (x1 >= 0 && x1 <= 1)
		x = x1
	else if (x2 >= 0 && x2 <= 1)
		x = x2
	else
		stop ("no max in interval")

	if (has.cx)
		x = cx1 + dx * x
	x
}

.intt.integral = function (cx1, cx2, cy1, cy2, ct1, ct2)
{   has.cx = (!missing (cx1) && !missing (cx2) )
	if (has.cx)
    {   dx = cx2 - cx1
        ct1 = dx * ct1
        ct2 = dx * ct2
    }
    params = .intt.params (cy1, cy2, ct1, ct2)
	y = sum (1 / (1:4) * params)
	if (has.cx)
		y = dx * y
	y
}

.intt.integral.2 = function (cx1, cx2, cy1, cy2, ct1, ct2, x)
	x * .intt.integral (cx1, cx2, cy1, cy2, ct1, ct2)

.spline.eval = function (nc, cx, cy, ct, x)
{	if (x < cx [1] || x > cx [nc])
		stop ("x value outside control points")
	nI = sum (cx <= x)
    if (x == cx [nI])
		cy [nI]
    else
		.intt (cx [nI], cx [nI + 1], cy [nI], cy [nI + 1], ct [nI], ct [nI + 1], x)
}

.test.nc = function (nc)
{	if (nc < 3)
		stop ("nc needs to be >= 3")
}

.slopes = function (n, x, y, correction=TRUE)
{	s = diff (y) / diff (x)
	t = c (s [1], (y [3:n] - y [1:(n-2)]) / (x [3:n] - x [1:(n-2)]), s [n - 1])
	if (correction)
	{	for (i in 2:(n-1) )
		{	s1 = sign (s [i])
			s2 = sign (t [i])
			s3 = sign (t [i + 1])
			if (s1 == s2 && s1 == s3)
			{	if (s [i] == 0)
					t [i] = 0
				else if (t [i] / s [i] > 3)
					t [i] = 3 * s [i]
			}
			s1 = sign (s [i - 1])
			s2 = sign (t [i - 1])
			s3 = sign (t [i])
			if (s1 == s2 && s1 == s3)
			{	if (s [i - 1] == 0)
					t [i] = 0
				else if (t [i] / s [i - 1] > 3)
					t [i] = 3 * s [i - 1]
			}
		}
	}
	t
}