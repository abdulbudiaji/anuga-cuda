<html>
    <head>
        <style type="text/css">
  pre { cursor: default; }
        </style>
    <script type="text/javascript">

// list of aliases
_aliases = {
<<<ALIASES>>>
};

// hack to handle silly events behavior 
_insideElem = false;

// colorize aliased accesses
function colorize(elem,color)
{
  var aliases = _aliases[elem.id];
	if(aliases)
	{
		for (k in aliases)
		{ 
			document.getElementById(k).setAttribute('style', 'background-color: ' + (color ? color : aliases[k]));
		}
	}
}

// mouse over event
function do_onmouseover(elem)
{
    if(!_insideElem)
    {
        colorize(elem,undefined);
        _insideElem = true;
    }
}

// moue leave event
function do_onmouseout(elem)
{
    _insideElem = false;
    colorize(elem,"transparent");   
}

    </script>
    </head>
    <body>
<!-- dpil code unparsed -->
        <pre>
<<<CODE>>>
        </pre>
    </body>
</html>
