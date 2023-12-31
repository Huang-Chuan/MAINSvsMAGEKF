function out1 = dbtheta_dpsi(in1,in2,in3)
%DBTHETA_DPSI
%    OUT1 = DBTHETA_DPSI(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    18-Oct-2023 09:31:39

% Calculates d(B(psi)*theta)/dPsi
%
%Inputs: dp d Ts omega theta
%
%
d1 = in2(1,:);
d2 = in2(2,:);
d3 = in2(3,:);
dpsi1 = in1(1,:);
dpsi2 = in1(2,:);
dpsi3 = in1(3,:);
dpsi4 = in1(4,:);
dpsi5 = in1(5,:);
dpsi6 = in1(6,:);
th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th9 = in3(9,:);
th10 = in3(10,:);
th11 = in3(11,:);
th12 = in3(12,:);
th13 = in3(13,:);
th14 = in3(14,:);
th15 = in3(15,:);
t2 = d1.*dpsi5;
t3 = d2.*dpsi4;
t4 = d1.*dpsi6;
t5 = d3.*dpsi4;
t6 = d2.*dpsi6;
t7 = d3.*dpsi5;
t8 = d1.*th4;
t9 = d2.*th6;
t10 = d3.*th7;
t11 = dpsi4.*th4;
t12 = dpsi5.*th6;
t13 = dpsi6.*th7;
t14 = d1.*2.0;
t15 = d2.*2.0;
t16 = d3.*2.0;
t17 = d1.*6.0;
t18 = d3.*4.0;
t19 = d2.*6.0;
t20 = d3.*6.0;
t21 = d1.^2;
t22 = d2.^2;
t23 = d3.^2;
t24 = dpsi1.*2.0;
t25 = dpsi2.*2.0;
t26 = dpsi3.*2.0;
t27 = dpsi1.*6.0;
t28 = dpsi2.*6.0;
t29 = dpsi3.*6.0;
t30 = dpsi1.^2;
t31 = dpsi2.^2;
t32 = dpsi3.^2;
t33 = dpsi4.^2;
t34 = dpsi5.^2;
t35 = dpsi6.^2;
t36 = th5.*2.0;
t37 = th8.*2.0;
t67 = d3.*dpsi6.*-2.0;
t71 = d2.*th8.*-2.0;
t72 = dpsi4.*th5.*-2.0;
t38 = t2.*2.0;
t39 = t3.*2.0;
t40 = t4.*2.0;
t41 = t5.*2.0;
t42 = t2.*4.0;
t43 = t3.*4.0;
t44 = t6.*2.0;
t45 = t7.*2.0;
t46 = dpsi6.*t16;
t47 = t2.*6.0;
t48 = t3.*6.0;
t49 = t4.*6.0;
t50 = t5.*6.0;
t51 = t6.*6.0;
t52 = t7.*6.0;
t53 = t14.*th5;
t54 = t15.*th8;
t55 = dpsi4.*t36;
t56 = dpsi5.*t37;
t57 = dpsi4.*t8;
t58 = dpsi5.*t9;
t59 = dpsi6.*t10;
t60 = -t2;
t62 = -t5;
t63 = -t6;
t61 = -t38;
t64 = -t41;
t65 = -t42;
t66 = -t44;
t68 = -t47;
t69 = -t50;
t70 = -t51;
t73 = dpsi4.*t53;
t74 = dpsi5.*t54;
t75 = d1+dpsi1+t7+t63;
t76 = d2+dpsi2+t4+t62;
t77 = d3+dpsi3+t3+t60;
t78 = d1.*t75;
t79 = d2.*t76;
t80 = d3.*t77;
t81 = dpsi6.*t77;
t82 = t75.^2;
t83 = t76.^2;
t84 = t77.^2;
t85 = t14.*t76;
t86 = t14.*t77;
t87 = t15.*t75;
t88 = t15.*t76;
t89 = t15.*t77;
t90 = t16.*t75;
t91 = t16.*t76;
t92 = t16.*t77;
t93 = t17.*t77;
t94 = t19.*t77;
t95 = t20.*t75;
t96 = t20.*t76;
t97 = dpsi4.*t77.*2.0;
t98 = dpsi5.*t77.*2.0;
t99 = dpsi6.*t76.*2.0;
t101 = dpsi4.*t77.*6.0;
t102 = dpsi5.*t77.*6.0;
t108 = t38.*t76;
t109 = t38.*t77;
t110 = t39.*t75;
t111 = t39.*t77;
t112 = t40.*t77;
t114 = t41.*t77;
t115 = t44.*t76;
t116 = t44.*t77;
t117 = t45.*t77;
t118 = t46.*t76;
t119 = t46.*t77;
t120 = t50.*t77;
t121 = t52.*t77;
t122 = t14+t24+t45+t66;
t123 = t15+t25+t40+t64;
t124 = t16+t26+t39+t61;
t125 = t17+t27+t52+t70;
t126 = t19+t28+t49+t69;
t127 = t20+t29+t48+t68;
t130 = d1.*dpsi4.*t77.*-2.0;
t131 = t5.*t77.*-2.0;
t132 = t6.*t76.*-2.0;
t133 = t67.*t77;
t160 = t75.*t76;
t100 = t81.*2.0;
t103 = dpsi6.*t80;
t104 = t84.*3.0;
t105 = -t97;
t106 = -t81;
t107 = dpsi4.*t86;
t113 = dpsi5.*t89;
t128 = -t82;
t129 = -t84;
t134 = d1.*t122;
t135 = d3.*t122;
t136 = dpsi4.*t122;
t137 = dpsi4.*t123;
t138 = dpsi4.*t124;
t139 = dpsi5.*t122;
t140 = dpsi5.*t123;
t141 = dpsi5.*t124;
t142 = dpsi6.*t122;
t143 = dpsi6.*t123;
t144 = dpsi6.*t124;
t145 = dpsi4.*t126;
t146 = dpsi5.*t125;
t147 = dpsi6.*t127;
t149 = t2.*t122;
t150 = t4.*t122;
t159 = t60.*t122;
t161 = t76.*t122;
t162 = t77.*t122;
t163 = t77.*t123;
t164 = t77.*t125;
t165 = t77.*t126;
t166 = t86+t90;
t167 = t89+t91;
t168 = t98+t122;
t148 = dpsi4.*t134;
t151 = dpsi6.*t135;
t152 = -t134;
t153 = -t136;
t154 = -t139;
t155 = -t140;
t156 = -t144;
t157 = -t146;
t158 = -t147;
t169 = t105+t123;
t170 = t124+t137;
t171 = t123+t142;
t172 = t127+t145;
t173 = t172.*th10;
t174 = t124+t154;
t175 = t127+t157;
t176 = -t173;
t177 = t175.*th15;
t178 = -t177;
et1 = -t11+t12-t36-t37-d2.*th10.*6.0-d3.*th9.*2.0-d1.*th12.*2.0-d1.*th15.*6.0-d2.*th14.*2.0-d3.*th13.*2.0-dpsi2.*th10.*6.0-dpsi3.*th9.*2.0-dpsi1.*th12.*2.0-dpsi1.*th15.*6.0-dpsi2.*th14.*2.0-dpsi3.*th13.*2.0-t3.*th9.*4.0-t4.*th10.*6.0+t5.*th10.*1.2e+1-t3.*th13.*2.0-t4.*th14.*2.0+t5.*th14.*4.0-t7.*th12.*4.0-t7.*th15.*1.2e+1+t38.*th9+t42.*th13+t44.*th12+t51.*th15-d1.*dpsi4.*th11+d2.*dpsi5.*th11-dpsi2.*dpsi4.*th9.*2.0-dpsi1.*dpsi4.*th11+dpsi2.*dpsi5.*th11-dpsi3.*dpsi5.*th12.*2.0-dpsi3.*dpsi5.*th15.*6.0-dpsi4.*t2.*th10.*6.0-dpsi4.*t4.*th9.*2.0+dpsi6.*t2.*th11-dpsi4.*t2.*th14.*2.0-dpsi5.*t3.*th12.*2.0+dpsi6.*t3.*th11;
et2 = dpsi5.*t5.*th11.*-2.0-dpsi5.*t3.*th15.*6.0-dpsi5.*t6.*th13.*2.0+dpsi5.*t24.*th13+dpsi4.*t29.*th10+dpsi4.*t26.*th14+dpsi4.*t41.*th9+dpsi5.*t38.*th12+dpsi4.*t39.*th14+dpsi4.*t48.*th10+dpsi5.*t45.*th13+dpsi5.*t47.*th15;
et3 = t58+t71-th2-d2.*th5.*4.0-d3.*th4.*2.0-d1.*th7-dpsi2.*th5.*2.0-dpsi3.*th4-dpsi1.*th7+t2.*th4-t3.*th4.*2.0-t4.*th5.*2.0+t5.*th5.*4.0+t6.*th7-t7.*th7.*2.0-t22.*th10.*9.0+t23.*th10.*9.0-t21.*th14-t22.*th14.*2.0+t23.*th14.*3.0-t31.*th10.*3.0+t32.*th10.*3.0-t30.*th14+t32.*th14+t2.^2.*th10.*3.0+t3.^2.*th10.*9.0-t4.^2.*th10.*3.0-t5.^2.*th10.*9.0+t2.^2.*th14+t3.^2.*th14.*3.0-t7.^2.*th14.*3.0-d2.*d3.*th9.*6.0-d1.*d2.*th12.*4.0-d1.*d3.*th11.*2.0-d1.*d2.*th15.*6.0-d2.*d3.*th13.*2.0-d2.*dpsi2.*th10.*1.2e+1-d2.*dpsi3.*th9.*4.0-d3.*dpsi2.*th9.*4.0;
et4 = d1.*dpsi2.*th12.*-2.0-d1.*dpsi3.*th11-d2.*dpsi1.*th12.*4.0-d3.*dpsi1.*th11.*2.0-d1.*dpsi1.*th14.*2.0+d3.*dpsi3.*th10.*1.2e+1-d2.*dpsi1.*th15.*6.0-d2.*dpsi2.*th14.*2.0-d2.*dpsi3.*th13.*2.0-dpsi2.*dpsi3.*th9.*2.0-dpsi1.*dpsi2.*th12.*2.0-dpsi1.*dpsi3.*th11+d1.*t2.*th11-d2.*t3.*th9.*6.0-d1.*t3.*th11.*2.0-d3.*t2.*th10.*1.2e+1-d2.*t4.*th10.*1.2e+1+d3.*t3.*th10.*3.6e+1-d3.*t4.*th9.*4.0-d1.*t4.*th12.*2.0+d1.*t5.*th12.*4.0-d2.*t3.*th13.*2.0-d3.*t2.*th14.*8.0+d2.*t6.*th12.*4.0+d3.*t3.*th14.*8.0-d2.*t7.*th12.*8.0-d3.*t7.*th11.*3.0-d2.*t7.*th15.*1.2e+1+d2.*t42.*th9+d2.*t42.*th13+dpsi1.*t2.*th11-dpsi2.*t3.*th9.*4.0-dpsi1.*t3.*th11.*2.0-dpsi3.*t2.*th10.*6.0;
et5 = dpsi2.*t4.*th10.*-6.0+dpsi3.*t3.*th10.*1.2e+1-dpsi3.*t4.*th9.*2.0-dpsi1.*t4.*th12.*2.0+dpsi2.*t5.*th10.*1.2e+1+dpsi3.*t5.*th9.*4.0+dpsi1.*t5.*th12.*4.0-dpsi3.*t2.*th14.*2.0+dpsi3.*t6.*th11-dpsi2.*t7.*th12.*4.0-dpsi3.*t7.*th11.*2.0-dpsi1.*t7.*th14.*4.0+dpsi3.*t18.*th14+dpsi5.*t22.*th11+dpsi3.*t43.*th14-t2.*t3.*th10.*1.2e+1-t2.*t5.*th9.*4.0-t3.*t4.*th9.*4.0-t2.*t3.*th14.*4.0+t4.*t5.*th10.*1.2e+1-t3.*t7.*th11.*4.0+t6.*t7.*th14.*4.0+t6.*t16.*th11+t5.*t20.*th9+t2.*t25.*th9+t6.*t19.*th15+t6.*t25.*th12+t6.*t24.*th14+t4.*t38.*th9+t6.*t39.*th11+t7.*t38.*th11+t6.*t40.*th12+t5.*t48.*th9+t7.*t50.*th12+t6.*t63.*th14+d2.*dpsi2.*dpsi5.*th11-d2.*dpsi3.*dpsi5.*th12.*2.0-d2.*dpsi3.*dpsi5.*th15.*6.0-d2.*dpsi5.*t3.*th12.*2.0;
et6 = d3.*dpsi6.*t2.*th12.*-4.0-d3.*dpsi6.*t3.*th12.*4.0-d2.*dpsi5.*t3.*th15.*6.0-d2.*dpsi5.*t6.*th13.*2.0+dpsi1.*dpsi5.*t15.*th13+dpsi5.*t2.*t15.*th12+dpsi5.*t7.*t15.*th13+dpsi5.*t2.*t19.*th15;
et7 = t53+t57+th3+d1.*th8.*4.0+d2.*th7+dpsi2.*th7+dpsi3.*th6-t2.*th6.*2.0+t3.*th6+t4.*th7-t5.*th7.*2.0-t6.*th8.*2.0+t7.*th8.*4.0+t16.*th6+t24.*th8+t22.*th12-t23.*th12.*3.0+t21.*th15.*9.0-t23.*th15.*9.0+t31.*th12-t32.*th12+t30.*th15.*3.0-t32.*th15.*3.0-t2.^2.*th12.*3.0-t3.^2.*th12+t4.^2.*th12-t2.^2.*th15.*9.0+t5.^2.*th12.*3.0-t3.^2.*th15.*3.0+t6.^2.*th15.*3.0+t7.^2.*th15.*9.0+d1.*d2.*th14.*4.0+d2.*dpsi3.*th11+d1.*dpsi1.*th15.*1.2e+1+d1.*dpsi2.*th14.*4.0+d1.*dpsi3.*th13.*4.0-d3.*dpsi3.*th12.*4.0-d3.*dpsi3.*th15.*1.2e+1+dpsi2.*dpsi3.*th11-d1.*t2.*th9.*2.0-d2.*t2.*th11.*2.0-d1.*t2.*th13.*6.0;
et8 = d1.*t5.*th10.*-1.2e+1+d2.*t3.*th11+d3.*t2.*th12.*8.0-d3.*t3.*th12.*8.0+d1.*t4.*th14.*4.0-d3.*t5.*th11.*3.0-d1.*t5.*th14.*8.0+d3.*t2.*th15.*3.6e+1-d2.*t4.*th15.*1.2e+1-d3.*t3.*th15.*1.2e+1-d2.*t6.*th14.*2.0-d3.*t6.*th13.*4.0+d2.*t7.*th14.*4.0+d3.*t14.*th9+d1.*t14.*th12+d2.*t17.*th10+d3.*t15.*th11+d3.*t17.*th13+d1.*t43.*th9+d1.*t43.*th13-dpsi2.*t2.*th11.*2.0-dpsi1.*t2.*th13.*4.0+dpsi2.*t3.*th11-dpsi3.*t3.*th12.*2.0+dpsi3.*t4.*th11-dpsi2.*t5.*th12.*4.0-dpsi3.*t5.*th11.*2.0-dpsi1.*t5.*th14.*4.0+dpsi3.*t2.*th15.*1.2e+1-dpsi3.*t3.*th15.*6.0-dpsi1.*t6.*th15.*6.0-dpsi2.*t6.*th14.*2.0-dpsi3.*t6.*th13.*2.0+dpsi1.*t7.*th15.*1.2e+1+dpsi2.*t7.*th14.*4.0+dpsi3.*t7.*th13.*4.0;
et9 = dpsi3.*t14.*th9+dpsi1.*t14.*th12+dpsi2.*t15.*th12+dpsi2.*t16.*th11+dpsi2.*t17.*th10+dpsi1.*t15.*th14+dpsi1.*t18.*th13+dpsi4.*t21.*th11+dpsi2.*t24.*th14+dpsi3.*t24.*th13+dpsi3.*t42.*th12-t2.*t4.*th11.*2.0-t3.*t5.*th11.*2.0+t2.*t3.*th15.*1.2e+1-t4.*t5.*th12.*4.0-t2.*t7.*th13.*6.0-t3.*t6.*th13.*2.0-t4.*t6.*th14.*2.0-t5.*t7.*th14.*6.0-t6.*t7.*th15.*1.2e+1+t4.*t16.*th11+t4.*t17.*th10+t3.*t24.*th13+t7.*t20.*th13+t4.*t25.*th12+t4.*t24.*th14+t3.*t42.*th12+t5.*t42.*th11+t6.*t42.*th13+t7.*t43.*th13+d1.*dpsi1.*dpsi4.*th11-d1.*dpsi3.*dpsi4.*th10.*6.0-d1.*dpsi3.*dpsi4.*th14.*2.0-d1.*dpsi4.*t3.*th10.*6.0-d1.*dpsi4.*t5.*th9.*2.0-d1.*dpsi4.*t3.*th14.*2.0+dpsi2.*dpsi4.*t14.*th9+dpsi4.*t4.*t14.*th9+dpsi4.*t2.*t17.*th10+dpsi4.*t2.*t14.*th14;
et10 = dpsi6.*t2.*t18.*th14+dpsi6.*t3.*t18.*th14;
mt1 = [-t12+t13+t37+th12.*(t98+t99)+th15.*(t102+t125)+t171.*th14+t174.*th13+th11.*(t81-dpsi5.*t76),th7-th13.*(t100+t153)+dpsi4.*th6-dpsi6.*th8.*2.0+t169.*th12-th14.*(t99-t122)+th11.*(t77+dpsi4.*t76)-th15.*(t101+dpsi6.*t125),t56+t178+th6-dpsi4.*th7+t168.*th13+th11.*(t76-dpsi4.*t77)-th12.*(t124+dpsi4.*t76.*2.0)-th14.*(t136-dpsi5.*t76.*2.0),th7+th9.*(t100+t155)+dpsi6.*t36-dpsi5.*th4+t168.*th14+t171.*th12+th11.*(t77-dpsi5.*t75)+th10.*(t102+dpsi6.*t126)];
mt2 = [t11-t13+t36-th14.*(t97+t142)+t170.*th9-th10.*(t101-t126)+th12.*(t122-t143)-th11.*(t81-dpsi4.*t75),t72+t176+th4-th12.*(t136+t155)+dpsi5.*th7+t169.*th9-t174.*th14+th11.*(t75+dpsi5.*t77),t56+t178+th6+th13.*(t122+t141)+th9.*(t141+t143)+th14.*(t140+t156)+dpsi5.*t36+dpsi6.*th4-t174.*th12+th11.*(t76+dpsi6.*t75)-th10.*(t147-dpsi5.*t126),t72+t176+th4-th13.*(t138+t142)-th12.*(t136+t156)-dpsi4.*th8.*2.0-dpsi6.*th6-t170.*th14+th9.*(t123-t138)+th11.*(t75-dpsi6.*t76)+th15.*(t147-dpsi4.*t125),et1+et2];
mt3 = [t9-t10+t74+th11.*(t79-t80+t6.*t75+t7.*t75)+th4.*(t6+t7)+th9.*(t133+dpsi5.*t167+t6.*t123)-th10.*(t121+dpsi6.*(t94+t96)-d2.*dpsi5.*t126)-th12.*(t151+t167+d2.*t154)-th14.*(t116+t117+t135+d2.*t155)+th13.*(t113+d2.*t122)-th15.*(t94+d2.*t157)-th5.*(t46-dpsi5.*t15)];
mt4 = [t59+th1-th5.*(t18+t26+t43+t61)-th8.*(t16+t26+t43+t61)-th10.*(t94+t96+t165-t5.*t77.*6.0+t3.*t126)-th15.*(t164-t6.*t77.*6.0+t3.*t125)-th9.*(-t83+t84+t92-d2.*t123+dpsi4.*t167)-th12.*(t135+t162-dpsi6.*t167+t3.*t122)-th14.*(t89+t131-t151+t163+t3.*t123)-th13.*(t84+t111+t128+t6.*t122)+th11.*(t103+t160+d2.*t75+t62.*t75+t63.*t76)+th6.*(d1+dpsi1+t7+t66)+th4.*(dpsi2+t4+t15+t64),et3+et4+et5+et6];
mt5 = [t59-th1+th5.*(t16+t26+t39+t65)+th8.*(t18+t26+t39+t65)+th10.*(t165+t49.*t77+t60.*t126)-th11.*(-t103+t160+d1.*t76+t4.*t75+t7.*t76)+th15.*(t93+t95+t121+t164+t60.*t125)+th12.*(t86+t117+t118+t159+t162)-th13.*(t80.*-2.0+t82+t129+t134+dpsi5.*t166)+th14.*(t91+t163+dpsi6.*t166+t60.*t123)-th9.*(t83+t109+t129+t4.*t123)-th4.*(d2+dpsi2+t40+t62)-th6.*(dpsi1+t14+t45+t63)];
mt6 = [-t8+t10+t73+th11.*(-t78+t80+t4.*t76+t5.*t76)+th6.*(t4+t5)+th14.*(t166+d1.*t137+t67.*t76)+th12.*(t91-t112+t148+t64.*t77)-th15.*(t120+dpsi6.*(t93+t95)-d1.*dpsi4.*t125)+th13.*(t133+t150+dpsi4.*t166)+th10.*(t93+d1.*t145)-th9.*(t130+d1.*t123)-th8.*(t46-dpsi4.*t14),et7+et8+et9+et10];
mt7 = [t58+t71+th2-th13.*(t89-d2.*dpsi5.*t75.*2.0)+th5.*(t4.*4.0+t15+t25+t64)+th9.*(t112+t163-t2.*t76.*2.0)+th10.*(t83.*3.0-t104+t47.*t77+t49.*t76)+th11.*(d1.*t77+dpsi5.*t79+t60.*t75+t63.*t77+t75.*t77)+th12.*(t85+t132+t150+t161-d2.*dpsi5.*t77.*2.0)-th15.*(dpsi5.*t94+t19.*t75)+th4.*(d3+dpsi3+t3+t61)+th14.*(t79.*-2.0+t82+t129+t134+t2.*t77.*2.0+t66.*t75)+th7.*(dpsi1+t7+t14+t66)];
mt8 = [t53+t57-th3+th15.*(t82.*-3.0+t104+t48.*t77+t51.*t75)+th10.*(t17.*t76-d1.*dpsi4.*t77.*6.0)-th8.*(t6.*-4.0+t14+t24+t45)-th13.*(t110+t162-t6.*t77.*2.0)-th14.*(t87+t107+t132+t150+t161)-th11.*(d2.*t77-dpsi4.*t78+t3.*t76+t4.*t77+t76.*t77)+th9.*(t86+dpsi4.*t85)-th6.*(d3+dpsi3+t39+t60)-th12.*(t79.*2.0+t83+t129+t152-t3.*t77.*2.0+t40.*t76)-th7.*(dpsi2+t15+t40+t62)];
mt9 = [t8-t9+th11.*(t78-t79+t2.*t77+t3.*t77)+th7.*(t2+t3)-th13.*(t87+t113)+th9.*(t85+t130)-th14.*(t86-t110+t159+dpsi5.*t79.*2.0)+th12.*(t89+t108-t148+t39.*t76)+d1.*t72+dpsi5.*t71+th15.*(t94-d2.*dpsi5.*t75.*6.0)-th10.*(t93+dpsi4.*t17.*t76)];
out1 = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9],3,6);
