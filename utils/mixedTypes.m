%To add custom types add a case statement with custom data type name. To
%invoke a specific type within another function use 
%T = mytypes(<data type name>) then cast as follows: x = cast(x0, 'like', T.x)
%To optimize the data type for a specific application use 'pg_solv_mex'.

function T = mixedTypes(dt) 
  switch(dt) 
    case 'double' 
      T.n = double([]);
      T.param = double([]);
      T.y = double([]);
      T.x = double([]);
    case 'single' 
      T.n = single([]); 
      T.param = single([]);
      T.y = single([]);
      T.x = single([]);
%     case 'fixed8' 
%       T.n = fi([],true,8,0); 
%       T.param = fi([],true,8,2);
%       T.x = fi([],true,8,4); 
%       T.y = fi([],true,8,4); 
    case 'fixed8' 
      T.n = fi([],true,8,0); 
      T.param = fi([],true,8,2);
      T.x = fi([],true,8,2); 
      T.y = fi([],true,8,2); 
    case 'fixed16' 
      T.n = fi([],true,16,0); 
      T.param = fi([],true,16,4);
      T.y = fi([],true,16,4);
      T.x = fi([],true,16,4);
    case 'fixed12' 
      T.n = fi([],true,12,2); 
      T.param = fi([],true,12,2);
      T.y = fi([],true,12,2);
      T.x = fi([],true,12,2);
    case 'fixed2' 
      T.n = fi([],true,2,0); 
      T.param = fi([],true,2,0);
      T.y = fi([],true,2,0);
      T.x = fi([],true,2,0);
    case 'fixed16override' 
      %use scaled doubles to detect potential overflows (see pg_solv_mex.m)
      T.n = fi([],true,16,0,'DataType', 'ScaledDouble');
      T.param = fi([],true,16,8,'DataType', 'ScaledDouble');      
      T.x = fi([],true,16,14,'DataType', 'ScaledDouble');
      T.y = fi([],true,16,10,'DataType', 'ScaledDouble');
     case 'fixed16+matchproc'
      %Types defined with fimath settings that match processor types 
      F = fimath(... 
      'RoundingMethod','Floor', ... 
      'OverflowAction','Wrap', ... 
      'ProductMode','KeepLSB', ... 
      'ProductWordLength',32, ... 
      'SumMode','KeepLSB', ... 
      'SumWordLength',32); 
      T.n = fi([],true,16,0,F);
      T.param = fi([],true,16,8,F);
      T.x = fi([],true,16,14,F);
      T.y = fi([],true,16,10,F);
  end 
end 