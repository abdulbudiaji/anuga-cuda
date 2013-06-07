require "caps"
require "bgf"

function bgf.escape(message)
  local tmp = string.gsub(message,"\n", "\\n")
  tmp = string.gsub(tmp,"\t", "     ")
  return tmp
end

function bgf._dumpComments(level, stmt)
   local indent = string.rep(" ",level*2)

   if (stmt:commentCount() > 0) then
     for j=0, stmt:commentCount()-1 do
       print (indent .. "  .cmnt \"" .. bgf.escape(stmt:comment(j)) .. "\"")
     end
   end
end

function bgf.dumpComments(stmt) 
   bgf._dumpComments(0,stmt)
end

function bgf._dumpStmt(level,stmt)

   local indent = string.rep(" ",level*2)

   repeat

      print(indent .. bgf.bf_type(stmt) .. " : (" .. stmt:id() .. ")")

      local name = stmt:name()
      if name then
         print(indent.."   -> name = "..name:identifier())
      end

      -- print comments
      bgf._dumpComments(level, stmt)

      local arg0 = stmt:operand0()
      local arg1 = stmt:operand1()
      local arg2 = stmt:operand2()

      bgf._dumpExpression(level+1, arg0)

      if (arg1) then
         bgf._dumpExpression(level+1, arg1)
      end

      if (arg2) then
         bgf._dumpExpression(level+1, arg2)
      end

      for j=0,stmt:childCount()-1 do
         local child = stmt:child(j)
         bgf._dumpStmt(level+1,child)
      end
      
      stmt = stmt:nextContainerInChain()
   until not stmt
end

function bgf.dumpStmt(stmt)
   bgf._dumpStmt(0,stmt)
end

function bgf._dumpExpression(level, expr)
   local indent = string.rep(" ",level*2)

   if (not expr) then
      print (indent.." | NULL")
      return
   end

   print(indent .. " | " .. bgf.bf_type(expr).. " : (" .. expr:id() .. ")    " .. bgf._dumpType(0, expr:type()) )

   local arg0 = expr:operand0()
   local arg1 = expr:operand1()

   bgf._dumpExpression(level+1, arg0)

   if (arg1) then
      bgf._dumpExpression(level+1, arg1)
   end
end

function bgf.dumpExpression(expr)
   bgf._dumpExpression(0,expr)
end

function bgf.dumpType(type)
   return bgf._dumpType(0,type)
end

function bgf._dumpType(level, type)
   local indent = string.rep(" ",level*2)
   
   if not type then return "" end

   local print_type = ""

   if (type:isReal()) then
     print_type = print_type .. "REAL"
   elseif type:isInteger() then
     print_type = print_type .. "INTEGER"     
   elseif (type:isLogical()) then
     print_type = print_type .. "LOGICAL"
   elseif (type:isCharacter()) then
     print_type = print_type .. "CHARACTER"
   elseif (type:isStructure()) then
     print_type = print_type .. "REAL" 
   end

   print_type = "meta=" .. print_type .. "(kind=" .. type:kind() .. ")"
   if (type:dimensionCount() > 0) then
     print_type = print_type .. "[" .. type:dimensionCount() .. "]"
   end

   return print_type

end


function bgf.dumpSymbols(symb)
   bgf._dumpSymbols(0, symb)
end

function bgf._dumpSymbols(level, symb)
   local indent = string.rep(" ",level*2)

   local print_symb = indent .. "SYMB #" .. symb:id() .. " : " .. symb:identifier() .. "    " .. symb:variant()
.. " : " .. bgf._dumpType(0,symb:metaType())
   
   if (symb:hasPARAMETER()) then print_symb = print_symb .. "  parameter" end
   if (symb:hasEXTERNAL()) then  print_symb = print_symb .. "  external" end

   print (print_symb)

   local child = symb:firstInScope()
   
   if child then
      repeat
        bgf._dumpSymbols(level+1, child)
        child = child:nextInScope()
      until not child
   end
end

function bgf.dump(file)
   for i=0,file:numberOfTopLevelStmt()-1 do
      local unit = file:getTopLevelStmt(i)
      bgf.dumpStmt(unit)
   end

   for i=0,file:numberOfTopLevelStmt()-1 do
      local unit = file:getTopLevelStmt(i)
      bgf.dumpSymbols(unit:name())
   end
end