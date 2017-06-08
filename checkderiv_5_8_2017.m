 close(figure(67));figure(67);
 hold on
 subplot(2,2,1)
 for pixels = what_bg_pxls;
     hold on
     plot(bgphaseuw_shft(:,pixels),'g');
 end
 plot(mean(bgphaseuw_shft(:,bg_ww),2),'k')
 grid on
%  xlim([580 4400]);
%  legend(['1','2','3','4','5'])
 clear pixels
 xlim([aaa bbb])
 
 subplot(2,2,2)
 for pixels = what_ref_pxls;
     hold on
     plot(refphaseuw_shft(:,pixels),'m');
 end
 plot(mean(refphaseuw_shft(:,ref_ww),2),'k')
 grid on
 clear pixels
 xlim([aaa bbb])
 
 subplot(2,2,3)
 for pixels = what_trans_pxls;
     hold on
     plot(transphaseuw_shft(:,pixels),'r');
 end
 plot(mean(transphaseuw_shft(:,trans_ww),2),'k')
 grid on
 clear pixels
 xlim([aaa bbb])
%  xlim([580 4400]);
 %%
 close(figure(57));figure(57);
 clear what_bg_pxls bg_diff_cnt chk_bg_diff_met chk_bg_diff
 clear what_ref_pxls ref_diff_cnt chk_ref_diff_met chk_ref_diff
 clear what_trans_pxls trans_diff_cnt chk_trans_diff_met chk_trans_diff
%  hold on

thrsh_bg=2.88;
thrsh_ref=2.75;%2.752;%2.96;
thrsh_trans=2.2;%2.416;
stopshft = 13300;%40000;
ofst = 1;%11500;

 subplot(3,1,1)
  bg_diff_cnt = 1;
 for pixels = 1:size(bg_pts,1);
% for pixels = 818:946;
    chk_bg_diff = sort(abs(diff(bgphaseuw_shft(shft+ofst:stopshft,pixels))),'descend');
    chk_bg_diff_met(pixels) = mean(chk_bg_diff(1:6));
    if chk_bg_diff_met(pixels) < thrsh_bg
        what_bg_pxls(bg_diff_cnt) = pixels;
        bg_diff_cnt = bg_diff_cnt+1;
    end
    clear chk_bg_diff
 end
 plot(chk_bg_diff_met,'go')
 hold on
 plot(what_bg_pxls,chk_bg_diff_met(what_bg_pxls),'gx')
 plot(what_bg_pxls,chk_bg_diff_met(what_bg_pxls),'ko')

 subplot(3,1,2)
 ref_diff_cnt = 1;
%  for pixels = 1:129;
 for pixels = 1:size(ref_pts,1);
    chk_ref_diff = sort(abs(diff(refphaseuw_shft(shft+ofst:stopshft,pixels))),'descend');
    chk_ref_diff_met(pixels) = mean(chk_ref_diff(1:6));
    if chk_ref_diff_met(pixels) < thrsh_ref
        what_ref_pxls(ref_diff_cnt) = pixels;
        ref_diff_cnt = ref_diff_cnt+1;
    end
    clear chk_ref_diff
 end
 plot(chk_ref_diff_met,'mo')
 hold on
 plot(what_ref_pxls,chk_ref_diff_met(what_ref_pxls),'mx')
 plot(what_ref_pxls,chk_ref_diff_met(what_ref_pxls),'ko')
 
 subplot(3,1,3)
 trans_diff_cnt = 1;
 for pixels = 1:size(trans_pts,1);
    chk_trans_diff = sort(diff(transphaseuw_shft(shft:stopshft,pixels)),'descend');
    chk_trans_diff_met(pixels) = mean(chk_trans_diff(1:6));
    if chk_trans_diff_met(pixels) < thrsh_trans
        what_trans_pxls(trans_diff_cnt) = pixels;
        trans_diff_cnt = trans_diff_cnt+1;
    end
    clear chk_trans_diff
 end
 plot(chk_trans_diff_met,'ko')
 hold on
 plot(what_trans_pxls,chk_trans_diff_met(what_trans_pxls),'rx')
 plot(what_trans_pxls,chk_trans_diff_met(what_trans_pxls),'ro')
% %  plot(mean(bgphaseuw_shft(:,1:323),2),'k')
% %  xlim([990 4400]);
% %  legend(['1','2','3','4','5'])