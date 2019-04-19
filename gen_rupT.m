% Coded by Song, S. (Feb. 2010)
% Computing rupture time dist. from local rupture velocity 
% assuming straight line between hypo and point of interest

function [trup, lrup] = gen_rupT(Vr,lx,lz)

nx = length(lx); 
nz = length(lz);

dx = lx(2) - lx(1); 
dz = lz(2) - lz(1);

[dummy1 hi] = min(abs(lx));
[dummy2 hk] = min(abs(lz));

% move hypo to the nearest node
lx = lx - lx(hi);
lz = lz - lz(hk);

for i=1:nx 
   for k=1:nz
      trup(k,i) = 0;
      lrup(k,i) = 0;

      if (i ~= hi) | (k ~= hk)
         azi = atan2(lz(k),lx(i));
              
         if (azi >= 0) & (azi < pi/2)

             ci = hi;
             ck = hk;

             cz = dz/2; 
             cx = dx/2;

             while ((ci <= i) & (ck <= k))
                cz_t = cz - tan(azi)*cx;
                if (cz_t > 0)
                   dl = cx/cos(azi);
                   lrup(k,i) = lrup(k,i) + dl;
                   trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                   ci = ci + 1;
                   cz = cz - tan(azi)*cx;
                   cx = dx;
                else
                   dl = cz/sin(azi);
                   lrup(k,i) = lrup(k,i) + dl;
                   trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                   ck = ck + 1;
                   cx = cx - cz/tan(azi);
                   cz = dz;
                end
            
                if (ci == i) & (ck == k)
                   if cx == dx
                      dl = cx/cos(azi)/2;
                      lrup(k,i) = lrup(k,i) + dl;
                      trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                   else
                      dl = cz/sin(azi)/2;
                      lrup(k,i) = lrup(k,i) + dl;
                      trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                   end    
                
                   break;
                end    
             end          

         elseif (azi >= pi/2) & (azi < pi)

            ci = hi;
            ck = hk;
            
            cz = dz/2; 
            cx = dx/2;

            azi = pi - azi;
            while ((ci >= i) & (ck <= k))
            
               cz_t = cz - tan(azi)*cx;
               if (cz_t > 0)
                  dl = cx/cos(azi);
                  lrup(k,i) = lrup(k,i) + dl;
                  trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  ci = ci - 1;
                  cz = cz - tan(azi)*cx;
                  cx = dx;
               else
                  dl = cz/sin(azi);
                  lrup(k,i) = lrup(k,i) + dl;
                  trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  ck = ck + 1;
                  cx = cx - cz/tan(azi);
                  cz = dz;
               end
            
               if (ci == i) & (ck == k)
                  if cx == dx
                     dl = cx/cos(azi)/2;
                     lrup(k,i) = lrup(k,i) + dl;
                     trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                   else
                     dl = cz/sin(azi)/2;
                     lrup(k,i) = lrup(k,i) + dl;
                     trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                   end    
                
                   break;
               end    
                       
            end
            
         elseif (azi < 0) & (azi >= -pi/2)

            ci = hi;
            ck = hk;
            
            cz = dz/2; 
            cx = dx/2;

            azi = abs(azi);
            while ((ci <= i) & (ck >= k))

               cz_t = cz - tan(azi)*cx;
               if (cz_t > 0)
                  dl = cx/cos(azi);
                  lrup(k,i) = lrup(k,i) + dl;
                  trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  ci = ci + 1;
                  cz = cz - tan(azi)*cx;
                  cx = dx;
               else
                  dl = cz/sin(azi);
                  lrup(k,i) = lrup(k,i) + dl;
                  trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  ck = ck - 1;
                  cx = cx - cz/tan(azi);
                  cz = dz;
               end

               if (ci == i) & (ck == k)
                  if cx == dx
                     dl = cx/cos(azi)/2;
                     lrup(k,i) = lrup(k,i) + dl;
                     trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  else
                     dl = cz/sin(azi)/2;
                     lrup(k,i) = lrup(k,i) + dl;
                     trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  end

                  break;
               end

            end

         else

            ci = hi;
            ck = hk;
            
            cz = dz/2; 
            cx = dx/2;

            azi = pi + azi;
            while ((ci >= i) & (ck >= k))

               cz_t = cz - tan(azi)*cx;
               if (cz_t > 0)
                  dl = cx/cos(azi);
                  lrup(k,i) = lrup(k,i) + dl;
                  trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  ci = ci - 1;
                  cz = cz - tan(azi)*cx;
                  cx = dx;
               else
                   dl = cz/sin(azi);
                   lrup(k,i) = lrup(k,i) + dl;
                   trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                   ck = ck - 1;
                   cx = cx - cz/tan(azi);
                   cz = dz;
               end

               if (ci == i) & (ck == k)
                  if cx == dx
                     dl = cx/cos(azi)/2;
                     lrup(k,i) = lrup(k,i) + dl;
                     trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  else
                     dl = cz/sin(azi)/2;
                     lrup(k,i) = lrup(k,i) + dl;
                     trup(k,i) = trup(k,i) + dl/Vr(ck,ci);
                  end

                  break;
               end

            end

         end  
      end
   end  % end of k-array
end  % end of i-array    


