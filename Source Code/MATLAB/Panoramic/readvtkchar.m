function [charbuff,whichterm]=readvtkchar(fid,asciiterm1,asciiterm2);
%
% reads chars from a vtk file until either of two asciiterm terminations are found
%
% asciiterm examples:
%   32: <SP> is a space, which is ASCII character 32
%   10: <NEW LINE> is a carriage return, which is ASCII character 10
%
% MWKay, 4/2005

havenothing=1;

while havenothing;
i=1;
gof=1;
whichterm=0;
while gof
  buff=fgets(fid,1);
  if (buff==asciiterm1) | (buff==asciiterm2) 
    gof=0;
    if (buff==asciiterm1)
      whichterm=1;
    else
      whichterm=2;
    end
  else
    charbuff(i)=buff;
    i=i+1;
    havenothing=0;
  end
end
end
