<html>
    <head>
        <style type="text/css">
  pre { cursor: default; }
        </style>
    <script type="text/javascript">

// list of defs
_defs = {
    <<<DEFS>>>
};

// list of uses
_uses = {
    <<<USES>>>
};

// hack to handle silly events behavior 
_insideElem = false;

// colorize defs and uses for the element
function colorize(elem,defColor,useColor)
{
    var uses = _defs[elem.id];
	if(uses)
	{
	    elem.setAttribute('style', 'background-color: ' + defColor);
		for (i=0;i<uses.length;i++)
		{
			document.getElementById(uses[i]).setAttribute('style', 'background-color: ' + useColor);
		}
	}

	var defs = _uses[elem.id];
	if(defs)
	{
	    elem.setAttribute('style', 'background-color: ' + useColor);
		for (i=0;i<defs.length;i++)
		{
			document.getElementById(defs[i]).setAttribute('style', 'background-color: ' + defColor);
		}
	}
}

// mouse over event
function do_onmouseover(elem)
{
    if(!_insideElem)
    {
        colorize(elem,"red","green");
        _insideElem = true;
    }
}

// moue leave event
function do_onmouseout(elem)
{
    _insideElem = false;
    colorize(elem,"transparent","transparent");   
}

    </script>
    </head>
    <body>
        <pre>
<!-- dpil code unparsed -->

<<<CODE>>>
        </pre>
    </body>
</html>
