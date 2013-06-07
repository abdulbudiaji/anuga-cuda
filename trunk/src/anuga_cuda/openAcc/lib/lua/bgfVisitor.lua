require "bgf"

-- Bf visitor class

function bgf.Visitor()
   local visitor = {}

   -- this is the hierachy class cache
   visitor._classnames_cache = {}

   -- this function is used to get class hierarchy
   visitor._getClassNames = 
      function(self,node)
         local id = node:classid()

         local result = self._classnames_cache[id]

         if result then
            return result
         end

         result = {}

         while id ~=0 do
            local classname = bgf.BfObject_classid2name(id)
            table.insert(result,classname)
            id = bgf.BfObject_classparent(id)
         end

         if node:isBfStatement() then
            table.insert(result,"BfStatement")
         elseif node:isBfExpression() then
            table.insert(result,"BfExpression")
         else
            error("Unexpected node in visitor")
         end

         -- visitor does not visit symbols and types
         -- if node:isBfSymbol() then
         --    table.insert(result,"BfSymbol")
         -- end

         -- if node:isBfType() then
         --    table.insert(result,"BfType")
         -- end

         self._classnames_cache[id] = result
         return result;
      end

   -- function used to call corresponding enterXXX
   visitor._doEnter =
      function(self,node)
         local names = self:_getClassNames(node);
         
         -- call the enter methods in the natural order
         for i=1,#names do
            local name = names[i]
            local enter = self.enter_table[name] 
            if enter then
               enter(self,node)
            end
         end

      end

   -- function used to call corresponding leaveXXX
   visitor._doLeave =
      function(self,node)
         local names = self:_getClassNames(node);

         -- call the leave methods in the reverse order
         for i=#names,1,-1 do
            local name = names[i]
            local leave = self.leave_table[name] 
            if leave then
               leave(self,node)
            end
         end

      end

   -- visit expr function
   visitor._visitExpr = 
      function(self, expr)

         self:_doEnter(expr)

         local arg0 = expr:operand0()
         local arg1 = expr:operand1()

         if arg0 then
            self:_visitExpr(arg0)
         end

         if arg1 then
            self:_visitExpr(arg1)
         end

         self:_doLeave(expr)

      end

   -- check if statment comments contains the 'str' string 
   visitor._checkComments =
      function(self,stmt,str)
         for i=0,stmt:commentCount()-1 do
            local comment = stmt:comment(i)
            if(string.match(comment,str)) then
               return true
            end
         end
         return false
      end

   -- visit stmt function
   visitor._visitStmt = 
      function(self,stmt,startLine,endLine)

         repeat           

            if self._partial then
               -- if visitor is finished return
               if self._finished == true then
                  return
               end

               -- start the visitor if the region comment is attached to the current statement 
               if not self._started then
                  if self:_checkComments(stmt,self._beginRegionString) then
                     self._started = true
                  end
               -- check for the of the region 
               else 
                  if self:_checkComments(stmt,self._endRegionString) then
                     self._finished = true
                     return
                  end
               end
            end
         
            if self._started then
               self:_doEnter(stmt)

               -- visit expr
               local arg0 = stmt:operand0()
               local arg1 = stmt:operand1()
               local arg2 = stmt:operand2()
            
               if arg0 then
                  self:_visitExpr(arg0)
               end

               if arg1 then
                  self:_visitExpr(arg1)
               end

               if arg2 then
                  self:_visitExpr(arg2)
               end
            end

            -- visit children
            for j=0,stmt:childCount()-1 do
               local child = stmt:child(j)
               self:_visitStmt(child)
            end

            if self._started then
               self:_doLeave(stmt)
            end

            -- go to next statement in chain
            stmt = stmt:nextContainerInChain()

         until not stmt
              
      end

   -- this function is the entry point of the visitor
   visitor.process = 
      function(self,file,regionId)

         -- build visit tables
         self.enter_table = {}
         self.leave_table = {}
        
         for k,v in pairs(self) do

            local enter = string.match(k,"enter(Bf.*)")
            if enter then
               if (not bgf.isBgfClassname(enter)) then
                  print("Warning: No class named '" .. enter .. "'")
               end
               self.enter_table[enter] = v
            end

            local leave = string.match(k,"leave(Bf.*)")

            if leave then
               if (not bgf.isBgfClassname(leave)) then
                  -- FIXME: try to display the user location of the leave function
                  print("Warning: No class named '" .. leave .. "'")
               end
               self.leave_table[leave] = v
            end
         end
        
         -- init search string in case we lookup a region
         if regionId then
            self._beginRegionString = "capstune begin region "..regionId 
            self._endRegionString = "capstune end region "..regionId
            self._started = false
            self._finished = false
            self._partial = true
         else
            self._beginRegionString = nil
            self._endRegionString = nil
            self._started = true
            self._finished = false
            self._partial = false
         end

         -- process top level statements
         for i=0,file:numberOfTopLevelStmt()-1 do
            
            if self._finished == true then
               return
            end
            
            local stmt = file:getTopLevelStmt(i)
            self:_visitStmt(stmt,startLine,endLine)
         end

      end

   return visitor
end

