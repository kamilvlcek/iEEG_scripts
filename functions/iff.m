function result = iff(test,trueVal,falseVal)
%IFF vyhodnoti test a vrati jeho hodnotu, true nebo false
%   z http://www.mathworks.com/matlabcentral/newsreader/view_thread/147044
         try
             if test
                 result = trueVal;
             else
                 result = falseVal;
             end
         catch
             result = false;
         end
 end



