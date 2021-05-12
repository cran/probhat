#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

"%$%" = function (object, name)
	attr (object, as.character (substitute (name) ) )
"%$%<-" = function (object, name, value)
	"attr<-" (object, as.character (substitute (name) ), value)

.LIST0 = function (..., .list.start.at=1)
{	sc = as.list (sys.call (-1) )[-(1:.list.start.at)]
	L = list (...)
	names.1 = names (sc)
	names.2 = as.character (sc)
	if (is.null (names.1) )
		names (L) = names.2
	else
	{	unnamed = (names.1 == "")
		names (L) = names.1
		names (L) [unnamed] = names.2 [unnamed]
	}
	L
}

.LIST = function (...)
	.LIST0 (...)

.EXTEND = function (object, class, ...)
{	cl = base::class (object)
	if (! missing (class) && ! is.null (class) )
		cl = c (class, cl)
	cl= list (cl)
	a1 = attributes (object)
	a1$class = NULL
	a2 = .LIST0 (..., .list.start.at=3)
	attributes (object) = c (class=cl, a1, a2)
	object
}

.UNPACK = function (x)
{	list2env (x, parent.frame (1) )
	invisible (NULL)
}

.THIS = function ()
	sys.function (-1)

.THAT = function ()
{	this = sys.function (-1)
	attributes (this)
}

.object.info = function (object, value, private.attributes, public.attributes, values, comment, n)
{	#print object's description
	c = class (object) [1]
	dims = .dim0 (object)
	cat (c, ", ", paste (dims, collapse=" * "), "\n", sep="")

	#print object's value
	if (value)
		.print.head (object, n)

	#print attributes
	. = attributes (object)
	com = .$comment
	if (!is.null (.) )
		. = .remove.special.attributes (.)
	if (!is.null (.) )
	{	n.attributes = length (.)
		names.attributes = names (.)
		if (n.attributes > 0)
		{	for (i in 1:n.attributes)
			{	char.1 = substring (names.attributes [i], 1, 1)
				if (private.attributes && char.1 == "." || public.attributes && char.1 != ".")
				{	x = . [[i]]
					c = class (x) [1]
					if (c == "matrix")
						c = paste (mode (x), c)
					dims = .dim0 (x)
					x.attributes = attributes (x)
					x.attributes = .remove.special.attributes (x.attributes)
					#print attribute description
					cat ("attr: ", names.attributes [i], ", ", c, ", ", paste (dims, collapse=" * "), sep="")
					if (!is.null (x.attributes) )
					{	names.x.attributes = names (x.attributes)
						char.1 = substring (names.x.attributes, 1, 1)
						n.x.private = sum (char.1 == ".")
						n.x.attributes = length (x.attributes)
						n.x.public = n.x.attributes - n.x.private
						n.print = 0
						if (private.attributes && public.attributes)
							n.print = n.x.attributes
						else if (private.attributes)
							n.print = n.x.private
						else if (public.attributes)
							n.print = n.x.public
						#print number of subattributes
						if (n.print > 0)
							cat (", (", n.print, ")", sep="")
					}
					cat ("\n")
					#print attribute value
					if (values)
						.print.head (. [[i]], n)
				}
			}
		}
	}
	#print comments
	if (comment && !is.null (com) )
		cat (com, sep="\n")
}

.object.summary = function (object, ...,
	value=TRUE, private.attributes=FALSE, public.attributes=TRUE, attribute.values=value, comment=TRUE,
	n=6)
	.object.info (object, value, private.attributes, public.attributes, attribute.values, comment, n)

.remove.special.attributes = function (.)
{	.$class = NULL
	.$comment = NULL
	.$srcref = NULL
	.$.Environment = NULL
	.$dim = NULL
	.$names = NULL
	.$dimnames = NULL
	.$row.names = NULL

	.
}

.dim0 = function (object)
{	dims = dim (object)
	if (is.null (dims) )
	{	if (is.function (object) )
		{	b = body (object)
			if (class (b) == "{")
				dims = length (format (body (object) ) )
			else
				dims = 1
		}
		else
			dims = length (object)
	}
	dims
}

.print.head = function (object, n)
{	if (isS4 (object) )
		print (object)
	else if (is.list (object) && !is.data.frame (object) )
		.print.head.list (object)
	else if (is.function (object) )
	{	attributes (object) = NULL
		h = head (object, n)
		for (l in h)
				cat (l, "\n", sep="")
	}
	else
	{	if (is.factor (object) )
			object = as.character (object)
		h = head (object, n)
		print (h, quote=FALSE)
	}
}

.print.head.list = function (object)
{	list.names = names (object)
	if (!is.null (list.names) )
	{	cat ("    %$% names\n")
		cat ("    (", paste (list.names, collapse=", "), ")\n", sep="")
	}
}
