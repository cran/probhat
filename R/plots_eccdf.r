plot.eccdf = function (x, with.points=FALSE, ...)
{   eccdf.f = x

	. = attributes (eccdf.f)
    cx = .$x
    cy = .$y
    x = seq (cx [1], cx [.$n], length.out=200)
	y = eccdf.f (x)
    plot (x, y, type="l", ...)
    if (with.points)
		points (cx, cy, pch=16)
}

lines.eccdf = function (x, ...)
{   eccdf.f = x

    . = attributes (eccdf.f)
    x = seq (.$x [1], .$x [.$n], length.out=200)
	y = eccdf.f (x)
	lines (x, y, type="l", ...)
}

plot.eccdf.inverse = function (x, with.points=FALSE, ...)
{   eccdf.f.inverse = x

    . = attributes (eccdf.f.inverse)
    y = seq (0, 1, length.out=200)
    x = eccdf.f.inverse (y)
    plot (y, x, type="l", ...)
    if (with.points)
        points (.$y, .$x, pch=16)
}

lines.eccdf.inverse = function (x, ...)
{   eccdf.f.inverse = x

    . = attributes (eccdf.f.inverse)
    y = seq (0, 1, length.out=200)
    x = eccdf.f.inverse (y)
    lines (y, x, type="l", ...)
}
