%% ===========================================================================
% SiStER_tangential_velocity_BCs
% Taper tangential velocity BC's so that they meet the velocity conditons
% on the adjoining boundaries gradually.  
%
% This was originally written for wedge models to avoid exporbitant pressure 
% build-up against the left side backstop do to the tangential motion on the base
% G.Ito 7/16
%% ===========================================================================


%-----------------------------------------------------------------------------
% If the BOTTOM is no slip
%-----------------------------------------------------------------------------
if (BC.bot(1)==0)
	if (length(BC.bot)<4); %if a tangential velocity is not defined then its zero
		BC.bot(4)=0;
	end

	BC.bot_profile=BC.bot(4)*ones(Nx,1); 

	if (BC.left(2)==0) %if left side is fixed
		jj=find(x<BC.bot_taperx(1));
		BC.bot_profile(jj)=x(jj).*(BC.bot(4)-BC.left(3))/(BC.bot_taperx(1));
	end

	if (BC.right(2)==0) %if right side is fixed
		jj=find(x>BC.bot_taperx(2));
		BC.bot_profile(jj)=BC.bot(4)+(x(jj)-BC.bot_taperx(2)).*(BC.right(3)-BC.bot(4))/(x(end)-BC.bot_taperx(2));
	end
end

%-----------------------------------------------------------------------------
% If the Bottom is no slip
%-----------------------------------------------------------------------------
if (BC.top(1)==0)
	if (length(BC.top)<4)
		BC.top(4)=0;
	end

	BC.top_profile=BC.top(4)*ones(Nx,1);

    if (BC.left(2)==0)  %if left side is fixed
    	jj=find(x<BC.top_taperx(1))
		BC.top_profile(jj)=BC.left(3) + x(jj).*(BC.top(4)-BC.left(3))/(BC.top_taperx(1));
	end
	if (BC.right(2)==0) %if right side is fixed
		jj=find(x>BC.top_taperx(2))
		BC.top_profile(jj)=BC.top(4) + (x(jj)-BC.top_taperx(2)).*(BC.right(3)-BC.top(4))/(x(end)-BC.top_taperx(2));
	end
end
