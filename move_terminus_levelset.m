function signeddistance_out = move_terminus_levelset(md, icelevelset_in, magnitude, direction, signeddistance_in_flag)

   % NOTE
   %  direction == -1: less ice (make signed distance values more positive)
   %            == +1: more ice (make signed distance values more negative)

   signeddistance_out = [];
   if size(icelevelset_in,1) ~= md.mesh.numberofvertices
      disp('ERROR: number of elements in levelset must equal number of vertices')
      return
   end

   % Initialize output
   signeddistance_out = nan * ones(size(icelevelset_in));

   % If a signed distance field is provided as input, use it directly
   if signeddistance_in_flag
      signeddistance_out = icelevelset_in;
   % If a mask is provided, convert to signed distances
   else
      % Create signed distance field
      disp('Converting levelsets to signed distance fields');
      for i=1:size(icelevelset_in,2)
         levelset = icelevelset_in(:,i);
         pos      = find(levelset<0);

         if exist('TEMP.exp','file'), delete('TEMP.exp'); end
         expcontourlevelzero(md, levelset, 0, 'TEMP.exp');
         % (TODO possibly) Check that there is only one zero-level contour
         levelset = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
         delete('TEMP.exp');
         levelset(pos) = -levelset(pos);
         signeddistance_out(:,i) = levelset;
      end
   end

   % Alter the signed distance field based on magnitude and direction
   switch sign(direction)
      case -1
         disp('Removing ice from the levelset(s)');
      case +1
         disp('Adding ice to the levelset(s)');
   end
   signeddistance_out = signeddistance_out - sign(direction) * magnitude;

end % main function

