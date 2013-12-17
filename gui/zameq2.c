#include "../zameq2.h"

#define MTR_URI "http://zamaudio.com/lv2/zameq2#"
#define MTR_GUI "ui"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "lv2/lv2plug.in/ns/extensions/ui/ui.h"

#define EQPOINTS 1000

#define LOGO_W (160.)
#define LOGO_H (30.)

#define PLOT_W (550.)
#define PLOT_H (380.)

typedef struct {
	LV2UI_Write_Function write;
	LV2UI_Controller controller;

	RobWidget *hbox, *ctable;

	RobTkLbl  *lbl_bw[4];
	RobTkSpin *knob_bw[4];
	RobTkLbl  *lbl_gain[4];
	RobTkSpin *knob_gain[4];
	RobTkLbl  *lbl_freq[4];
	RobTkSpin *knob_freq[4];
	RobTkSep  *sep[3];

	//RobTkDarea *darea;
	RobTkXYp  *xyp;
	RobTkImg  *logo;

	RobTkLbl  *lbl_ingain;
	RobTkSpin *knob_ingain;
	RobTkLbl  *lbl_outgain;
	RobTkSpin *knob_outgain;
	
	cairo_surface_t *frontface;
	cairo_surface_t *eqcurve;
	
	float eqx[EQPOINTS];
	float eqy[EQPOINTS];

	bool disable_signals;

} ZamEQ2_UI;

static inline double
to_dB(double g) {
	return (20.*log10(g));
}

static inline double
from_dB(double gdb) {
	return (exp(gdb/20.*log(10.)));
}

static inline double
sanitize_denormal(double value) {
	if (!isnormal(value)) {
		return (0.);
	}
	return value;
}

static void
peq(double G0, double G, double GB, double w0, double Dw,
        double *a0, double *a1, double *a2, double *b0, double *b1, double *b2, double *gn) {

	double F,G00,F00,num,den,G1,G01,G11,F01,F11,W2,Dww,C,D,B,A;
	F = fabs(G*G - GB*GB);
	G00 = fabs(G*G - G0*G0);
	F00 = fabs(GB*GB - G0*G0);
	num = G0*G0 * (w0*w0 - M_PI*M_PI)*(w0*w0 - M_PI*M_PI)
		+ G*G * F00 * M_PI*M_PI * Dw*Dw / F;
	den = (w0*w0 - M_PI*M_PI)*(w0*w0 - M_PI*M_PI)
		+ F00 * M_PI*M_PI * Dw*Dw / F;
	G1 = sqrt(num/den);
	G01 = fabs(G*G - G0*G1);
	G11 = fabs(G*G - G1*G1);
	F01 = fabs(GB*GB - G0*G1);
	F11 = fabs(GB*GB - G1*G1);
	W2 = sqrt(G11 / G00) * tan(w0/2.)*tan(w0/2.);
	Dww = (1.+ sqrt(F00 / F11) * W2) * tan(Dw/2.);
	C = F11 * Dww*Dww - 2. * W2 * (F01 - sqrt(F00 * F11));
	D = 2. * W2 * (G01 - sqrt(G00 * G11));
	A = sqrt((C + D) / F);
	B = sqrt((G*G * C + GB*GB * D) / F);
	*gn = G1;
	*b0 = (G1 + G0*W2 + B) / (1. + W2 + A);
	*b1 = -2.*(G1 - G0*W2) / (1. + W2 + A);
	*b2 = (G1 - B + G0*W2) / (1. + W2 + A);
	*a0 = 1.;
	*a1 = -2.*(1. - W2) / (1. + W2 + A);
	*a2 = (1. + W2 - A) / (1. + W2 + A); 

	*b1 = sanitize_denormal(*b1); 
	*b2 = sanitize_denormal(*b2);
	*a0 = sanitize_denormal(*a0);
	*a1 = sanitize_denormal(*a1);
	*a2 = sanitize_denormal(*a2);
	*gn = sanitize_denormal(*gn);
	if (!isnormal(*b0)) { *b0 = 1.; }
}

static void calceqcurve(float val[], float x[], float y[])
{
	for (uint32_t i = 0; i < EQPOINTS; ++i) {
		x[i] = i/(float)EQPOINTS;
		double L,M,N,O,P,Q,R;
		double complex H;
		double complex expiw = cos(-(i+0.0005)*M_PI/EQPOINTS*20./48.) + I*sin(-(i+0.0005)*M_PI/EQPOINTS*20./48.);
		double complex exp2iw = cos(-2*(i+0.0005)*M_PI/EQPOINTS*20./48.)+ I*sin(-2*(i+0.0005)*M_PI/EQPOINTS*20./48.);
		double freqH, phaseH;
		double dcgain = 1.f;

		double qq1 = pow(2.0, 1.0/val[4])/(pow(2.0, val[4]) - 1.0); //q from octave bw
		double boost1 = from_dB(val[8]);
		double fc1 = val[0] / 48000.;
		double w01 = fc1 * 2. * M_PI;
		double bwgain1 = sqrt(boost1)*sqrt(dcgain);
		double bw1 = fc1 / qq1;

		peq(1.0, boost1, bwgain1, w01, bw1, &P, &Q, &R, &M, &N, &O, &L);
		H = (M + N*expiw + O*exp2iw)/(P + Q*expiw + R*exp2iw);
		freqH = sqrt(creal(H)*creal(H)+cimag(H)*cimag(H));
		//phaseH = carg(H);

		y[i] = (to_dB(freqH)/70. + 5./7.) ;
		//printf("%.4f\n",y[i]);
	}
}

#include "gui/img/logo.c"

static void render_frontface(ZamEQ2_UI* ui) {

	cairo_t *cr;
	robtk_xydraw_set_surface(ui->xyp, NULL);
	ui->eqcurve = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, PLOT_W, PLOT_H);
	cr = cairo_create (ui->eqcurve);

	CairoSetSouerceRGBA(c_blk);
	cairo_set_operator (cr, CAIRO_OPERATOR_SOURCE);
	cairo_rectangle (cr, 0, 0, PLOT_W, PLOT_H);
	cairo_fill (cr);

	cairo_save(cr);
	rounded_rectangle (cr, 10, 10, PLOT_W - 20, PLOT_H - 20, 10);
	CairoSetSouerceRGBA(c_blk);
	cairo_fill_preserve(cr);
	cairo_clip(cr);
	cairo_set_operator (cr, CAIRO_OPERATOR_OVER);
	
	robtk_xydraw_set_surface(ui->xyp, ui->eqcurve);
	calceqcurve(ui->eqx, ui->eqx, ui->eqy);
	robtk_xydraw_set_points(ui->xyp, EQPOINTS, ui->eqx, ui->eqy);

}

static bool expose_event(RobWidget* handle, cairo_t* cr, cairo_rectangle_t *ev)
{
	ZamEQ2_UI* ui = (ZamEQ2_UI*) GET_HANDLE(handle);

	/* limit cairo-drawing to exposed area */
	cairo_rectangle (cr, ev->x, ev->y, ev->width, ev->height);
	cairo_clip(cr);
	cairo_set_source_surface(cr, ui->frontface, 0, 0);
	cairo_paint (cr);

	cairo_set_operator (cr, CAIRO_OPERATOR_OVER);

	return TRUE;
}

static void expose_plot(cairo_t* cr, void *smth)
{
}

static void xy_clip_fn(cairo_t *cr, void *data)
{
	ZamEQ2_UI* ui = (ZamEQ2_UI*) (data);
	rounded_rectangle(cr, 10, 10, PLOT_W-20, PLOT_H-20, 10);
	cairo_clip(cr);
	cairo_set_source_surface(cr, ui->eqcurve, 0, 0);
	cairo_paint(cr);
}

static bool cb_set_knobs (RobWidget* handle, void *data) {
	ZamEQ2_UI* ui = (ZamEQ2_UI*) (data);
	float val[12];
	if (ui->disable_signals) return TRUE;
	
	val[0] = robtk_spin_get_value(ui->knob_freq[0]);
	ui->write(ui->controller, ZAMEQ2_FREQ1, sizeof(float), 0, (const void*) &val[0]);
	val[1] = robtk_spin_get_value(ui->knob_freq[1]);
	ui->write(ui->controller, ZAMEQ2_FREQ2, sizeof(float), 0, (const void*) &val[1]);
	val[2] = robtk_spin_get_value(ui->knob_freq[2]);
	ui->write(ui->controller, ZAMEQ2_FREQL, sizeof(float), 0, (const void*) &val[2]);
	val[3] = robtk_spin_get_value(ui->knob_freq[3]);
	ui->write(ui->controller, ZAMEQ2_FREQH, sizeof(float), 0, (const void*) &val[3]);
	val[4] = robtk_spin_get_value(ui->knob_bw[0]);
	ui->write(ui->controller, ZAMEQ2_Q1, sizeof(float), 0, (const void*) &val[4]);
	val[5] = robtk_spin_get_value(ui->knob_bw[1]);
	ui->write(ui->controller, ZAMEQ2_Q2, sizeof(float), 0, (const void*) &val[5]);
	val[6] = robtk_spin_get_value(ui->knob_bw[2]);
	ui->write(ui->controller, ZAMEQ2_SLOPEDBL, sizeof(float), 0, (const void*) &val[6]);
	val[7] = robtk_spin_get_value(ui->knob_bw[3]);
	ui->write(ui->controller, ZAMEQ2_SLOPEDBH, sizeof(float), 0, (const void*) &val[7]);
	val[8] = robtk_spin_get_value(ui->knob_gain[0]);
	ui->write(ui->controller, ZAMEQ2_BOOSTDB1, sizeof(float), 0, (const void*) &val[8]);
	val[9] = robtk_spin_get_value(ui->knob_gain[1]);
	ui->write(ui->controller, ZAMEQ2_BOOSTDB2, sizeof(float), 0, (const void*) &val[9]);
	val[10] = robtk_spin_get_value(ui->knob_gain[2]);
	ui->write(ui->controller, ZAMEQ2_BOOSTDBL, sizeof(float), 0, (const void*) &val[10]);
	val[11] = robtk_spin_get_value(ui->knob_gain[3]);
	ui->write(ui->controller, ZAMEQ2_BOOSTDBH, sizeof(float), 0, (const void*) &val[11]);

	calceqcurve(val, ui->eqx, ui->eqy);
	robtk_xydraw_set_points(ui->xyp, EQPOINTS, ui->eqx, ui->eqy);

	return TRUE;
}

static void ui_disable(LV2UI_Handle handle)
{
}

static void ui_enable(LV2UI_Handle handle)
{
}

static RobWidget * toplevel(ZamEQ2_UI* ui, void * const top)
{

        ui->frontface = NULL;
	ui->hbox = rob_hbox_new(FALSE, 2);
        robwidget_make_toplevel(ui->hbox, top);
        ROBWIDGET_SETNAME(ui->hbox, "ZamEQ2");

        ui->ctable = rob_table_new(/*rows*/15, /*cols*/ 3, FALSE);
	ui->sep[0] = robtk_sep_new(TRUE);
	ui->sep[1] = robtk_sep_new(TRUE);
	ui->sep[2] = robtk_sep_new(TRUE);

	ui->knob_gain[0]   = robtk_spin_new(-50, 20, .1);
	ui->knob_gain[1]   = robtk_spin_new(-50, 20, .1);
	ui->knob_gain[2]   = robtk_spin_new(-50, 20, .1);
	ui->knob_gain[3]   = robtk_spin_new(-50, 20, .1);
	ui->knob_freq[0]   = robtk_spin_new(10, 20000, 1); // TODO log-map
	ui->knob_freq[1]   = robtk_spin_new(10, 20000, 1); // TODO log-map
	ui->knob_freq[2]   = robtk_spin_new(10, 20000, 1); // TODO log-map
	ui->knob_freq[3]   = robtk_spin_new(10, 20000, 1); // TODO log-map
	ui->knob_bw[0]   = robtk_spin_new(0.1, 6, .1);
	ui->knob_bw[1]   = robtk_spin_new(0.1, 6, .1);
	ui->knob_bw[2]   = robtk_spin_new(0.1, 1.1, .1);
	ui->knob_bw[3]   = robtk_spin_new(0.1, 1.1, .1);

	ui->lbl_gain[0] = robtk_lbl_new("              dB");
	ui->lbl_gain[1] = robtk_lbl_new("              dB");
	ui->lbl_gain[2] = robtk_lbl_new("              dB");
	ui->lbl_gain[3] = robtk_lbl_new("              dB");
	ui->lbl_freq[0] = robtk_lbl_new("LMF");
	ui->lbl_freq[1] = robtk_lbl_new("MF");
	ui->lbl_freq[2] = robtk_lbl_new("LF");
	ui->lbl_freq[3] = robtk_lbl_new("HF");
	ui->lbl_bw[0] = robtk_lbl_new("              Bw");
	ui->lbl_bw[1] = robtk_lbl_new("              Bw");
	ui->lbl_bw[2] = robtk_lbl_new("              Bw");
	ui->lbl_bw[3] = robtk_lbl_new("              Bw");

	//ui->darea = robtk_darea_new(PLOT_W,PLOT_H, &expose_plot, ui);

	ui->xyp = robtk_xydraw_new(PLOT_W, PLOT_H);
	robtk_xydraw_set_alignment(ui->xyp, 0, 0);
	robtk_xydraw_set_linewidth(ui->xyp, 2.5);
	robtk_xydraw_set_drawing_mode(ui->xyp, RobTkXY_yraw_line);
	robtk_xydraw_set_mapping(ui->xyp, 1., 0., 1., 0.);
	robtk_xydraw_set_area(ui->xyp, 10, 10, PLOT_W-20, PLOT_H-20);
	robtk_xydraw_set_clip_callback(ui->xyp, xy_clip_fn, ui);
	robtk_xydraw_set_color(ui->xyp, 1.0, .2, .0, 1.0);

	ui->logo = robtk_img_new(LOGO_W, LOGO_H, 4, img_logo.pixel_data);
	robtk_img_set_alignment(ui->logo, 0, 0);

	int row = 0;
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_bw[2]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_bw[2]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_freq[2]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_freq[2]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_gain[2]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_gain[2]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;

	rob_table_attach(ui->ctable, robtk_sep_widget(ui->sep[0]),
		0, 2, row, row+1, 2, 2, RTK_EXANDF, RTK_SHRINK);
	row++;

	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_bw[0]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_bw[0]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_freq[0]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_freq[0]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_gain[0]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_gain[0]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;

	rob_table_attach(ui->ctable, robtk_sep_widget(ui->sep[1]),
		0, 2, row, row+1, 2, 2, RTK_EXANDF, RTK_SHRINK);
	row++;

	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_bw[1]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_bw[1]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_freq[1]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_freq[1]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_gain[1]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_gain[1]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;

	rob_table_attach(ui->ctable, robtk_sep_widget(ui->sep[2]),
		0, 2, row, row+1, 2, 2, RTK_EXANDF, RTK_SHRINK);
	row++;

	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_bw[3]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_bw[3]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_freq[3]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_freq[3]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_lbl_widget(ui->lbl_gain[3]),
		0, 1, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	rob_table_attach(ui->ctable, robtk_spin_widget(ui->knob_gain[3]),
		1, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);
	row++;
	rob_table_attach(ui->ctable, robtk_img_widget(ui->logo),
		0, 2, row, row+1, 0, 0, RTK_EXANDF, RTK_SHRINK);

#define SPIN_DFTNVAL(SPB, VAL) \
	robtk_spin_set_default(SPB, VAL); \
	robtk_spin_set_value(SPB, VAL);

	SPIN_DFTNVAL(ui->knob_gain[0], 0)
	SPIN_DFTNVAL(ui->knob_gain[1], 0)
	SPIN_DFTNVAL(ui->knob_gain[2], 0)
	SPIN_DFTNVAL(ui->knob_gain[3], 0)
	SPIN_DFTNVAL(ui->knob_freq[0], 400)
	SPIN_DFTNVAL(ui->knob_freq[1], 1500)
	SPIN_DFTNVAL(ui->knob_freq[2], 250)
	SPIN_DFTNVAL(ui->knob_freq[3], 9000)
	SPIN_DFTNVAL(ui->knob_bw[0], 1)
	SPIN_DFTNVAL(ui->knob_bw[1], 1)
	SPIN_DFTNVAL(ui->knob_bw[2], 1)
	SPIN_DFTNVAL(ui->knob_bw[3], 1)

	robtk_spin_set_callback(ui->knob_freq[0], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_freq[1], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_freq[2], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_freq[3], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_bw[0], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_bw[1], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_bw[2], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_bw[3], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_gain[0], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_gain[1], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_gain[2], cb_set_knobs, ui);
	robtk_spin_set_callback(ui->knob_gain[3], cb_set_knobs, ui);

	rob_hbox_child_pack(ui->hbox, robtk_xydraw_widget(ui->xyp), FALSE, FALSE);
        rob_hbox_child_pack(ui->hbox, ui->ctable, FALSE, FALSE);
	
	return ui->hbox;
}

/******************************************************************************
 * LV2
 */

static LV2UI_Handle
instantiate(
		void* const               ui_toplevel,
		const LV2UI_Descriptor*   descriptor,
		const char*               plugin_uri,
		const char*               bundle_path,
		LV2UI_Write_Function      write_function,
		LV2UI_Controller          controller,
		RobWidget**               widget,
		const LV2_Feature* const* features)
{
	ZamEQ2_UI* ui = (ZamEQ2_UI*)calloc(1, sizeof(ZamEQ2_UI));

	if (!ui) {
		fprintf(stderr, "ZamEQ2 UI: out of memory\n");
		return NULL;
	}

	*widget = NULL;

	/* initialize private data structure */
	ui->write      = write_function;
	ui->controller = controller;

	*widget = toplevel(ui, ui_toplevel);
	render_frontface(ui);
	return ui;
}

static enum LVGLResize
plugin_scale_mode(LV2UI_Handle handle)
{
	return LVGL_LAYOUT_TO_FIT;
}

static void
cleanup(LV2UI_Handle handle)
{
	ZamEQ2_UI* ui = (ZamEQ2_UI*)handle;
	ui_disable(ui);

	robtk_xydraw_set_surface(ui->xyp, NULL);
	cairo_surface_destroy (ui->eqcurve);
	robtk_xydraw_destroy(ui->xyp);
	
	for (uint32_t i = 0; i < 4; ++i) {
		robtk_lbl_destroy(ui->lbl_gain[i]);
		robtk_lbl_destroy(ui->lbl_bw[i]);
		robtk_lbl_destroy(ui->lbl_freq[i]);
		robtk_spin_destroy(ui->knob_gain[i]);
		robtk_spin_destroy(ui->knob_bw[i]);
		robtk_spin_destroy(ui->knob_freq[i]);
	}

	//robtk_spin_destroy(ui->knob_ingain);
	//robtk_lbl_destroy(ui->lbl_ingain);
	//robtk_spin_destroy(ui->knob_outgain);
	//robtk_lbl_destroy(ui->lbl_outgain);
	robtk_img_destroy(ui->logo);

	rob_box_destroy(ui->hbox);
	cairo_surface_destroy(ui->frontface);

	free(ui);
}

static void
port_event(LV2UI_Handle handle,
		uint32_t     port_index,
		uint32_t     buffer_size,
		uint32_t     format,
		const void*  buffer)
{
	ZamEQ2_UI* ui = (ZamEQ2_UI*)handle;
	

	if (format != 0) return;
	const float v = *(float *)buffer;
	switch (port_index) {
		case ZAMEQ2_BOOSTDB1:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_gain[0], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_Q1:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_bw[0], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_FREQ1:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_freq[0], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_BOOSTDB2:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_gain[1], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_Q2:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_bw[1], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_FREQ2:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_freq[1], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_BOOSTDBL:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_gain[2], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_SLOPEDBL:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_bw[2], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_FREQL:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_freq[2], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_BOOSTDBH:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_gain[3], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_SLOPEDBH:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_bw[3], v);
			ui->disable_signals = false;
			break;
		case ZAMEQ2_FREQH:
			ui->disable_signals = true;
			robtk_spin_set_value(ui->knob_freq[3], v);
			ui->disable_signals = false;
			break;
		default:
			return;
	}
}

static const void*
extension_data(const char* uri)
{
	return NULL;
}

