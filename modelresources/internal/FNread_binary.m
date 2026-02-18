function V=read_binary(file,dims,prec,machineformat,iRec,recSize,numRecs);
% Generic function for reading binary files
% usage: V=read_binary(file,dims,prec,machineformat,iRec,recSize,numRecs);
% Defaults: 
%          prec = 'float64' (float64 = double precision, float32)
%          machineformat='ieee-be'; 
%          if number of arguments < 5, entire file is read
%          otherwise, the function reads numRecs record(s) of length 
%          recSize (of data type given by prec) starting at 
%          record iRec. If numRecs is omitted, only 1 record 
%          is read. If numRecs=Inf, all records (starting at iRec)
%          are read.
% Note: By default, the output array is of class 'double'. To force 
% the output array to have the SAME class (precision) as the input, 
% prepend a '*' before the precision argument, e.g., '*float32' or 
% '*float64' (the default, so not really necessary). More generally, 
% one can prescribe arbitrary precisions for the input data and 
% output variable by setting prec to 'source=>destination', e.g., 
% 'float32=>float64' or 'float64=>float32'.
% All this assumes that your version of MATLAB supports this.

% Samar Khatiwala (spk@ldeo.columbia.edu)

if nargin<2
  dims=[];
end  
if nargin<3
   prec='float64';
end
if nargin<4
   machineformat='ieee-be';
end
if isempty(prec)
  prec='float64';
end
if isempty(machineformat)
   machineformat='ieee-be';
end 

if nargin<7
  numRecs=1;
end

if nargin<5
  readSize=Inf;
  skipBytes=0;
else
  if ~isempty(findstr('64',prec))
    skipBytes=8*(iRec-1)*recSize;
  elseif ~isempty(findstr('32',prec)) 
    skipBytes=4*(iRec-1)*recSize;
  else    
    error('Unknown precision')
  end
  readSize=numRecs*recSize;
end
%skipBytes

fid=fopen(file,'r',machineformat);
status=fseek(fid,skipBytes,0);
V=fread(fid,readSize,prec); 
fclose(fid);

if ~isempty(dims)
  if length(dims)>1
	V=reshape(V,dims);
  end
end

