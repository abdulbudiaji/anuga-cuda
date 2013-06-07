PROGRAM=(<memory-spaces end ="\n\n"><statements end = "\n\n">)

DECLARE = <name>: <type>

TYPEDEF = type <name>: <type>

PARAMETERS = <statements separator=", ">

BASE_FUNCTION = <func_memspace end=" "><func_qualifier> <name>(<parameters>)<type begin=" : "><statements begin="\n">

INLINE_BLOCK = <statements separator = "; ">

BLOCK = ({\n<statements indent = "  " separator = "\n" epilogue = "\n">})

RETURN = return<value begin=" ">

PRAGMA = #<prefix><value>

WHILE = while <condition>\n<statements>

DOWHILE = do\n<statements>\nwhile <condition>

IF = if <condition>\n<statements separator = "\nelse\n">

FOR = for (<init> : <condition> : <cycle>)\n<statements>

LOOP = <kind end=" ">loop <var> : <type> in [<bounds separator=" : ">]\n<statements>

GOTO = goto <name>

LABEL = <name>::\n

BREAK = break

CONTINUE = continue

MEMSPACE = memspace <name> { <ptrsize> }

SWITCH = switch <condition>\n{\n<statements separator="\n">\n}

CASE_LIST = case <values separator=", ">
