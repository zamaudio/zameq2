#include "../zameq2.h"

#define MTR_URI "http://zamaudio.com/lv2/zameq2#"
#define MTR_GUI "ui"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lv2/lv2plug.in/ns/extensions/ui/ui.h"

/* widget, window size */
#define DAWIDTH  (620.)
#define DAHEIGHT (380.)
#define EQPOINTS 100

struct MyGimpImage {
        unsigned int   width;
        unsigned int   height;
        unsigned int   bytes_per_pixel;
        unsigned char  pixel_data[];
};

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

	RobTkXYp  *xyp;

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

/* load gimp-exported .c image into cairo surface */
static void img2surf (struct MyGimpImage const * img, cairo_surface_t **s, unsigned char **d) {
        unsigned int x,y;
        unsigned int stride = cairo_format_stride_for_width (CAIRO_FORMAT_ARGB32, img->width);

        (*d) = (unsigned char*) malloc (stride * img->height);
        (*s) = cairo_image_surface_create_for_data(*d,
                        CAIRO_FORMAT_ARGB32, img->width, img->height, stride);

        cairo_surface_flush (*s);
        for (y = 0; y < img->height; ++y) {
                const int y0 = y * stride;
                const int ys = y * img->width * img->bytes_per_pixel;
                for (x = 0; x < img->width; ++x) {
                        const int xs = x * img->bytes_per_pixel;
                        const int xd = x * 4;

                        if (img->bytes_per_pixel == 3) {
                        	(*d)[y0 + xd + 3] = 0xff;
                        } else {
                        	(*d)[y0 + xd + 3] = img->pixel_data[ys + xs + 3]; // A
                        }
                        (*d)[y0 + xd + 2] = img->pixel_data[ys + xs];     // R
                        (*d)[y0 + xd + 1] = img->pixel_data[ys + xs + 1]; // G
                        (*d)[y0 + xd + 0] = img->pixel_data[ys + xs + 2]; // B
                }
        }
        cairo_surface_mark_dirty (*s);
}

//#include "gui/img/eq2.c"

static void render_frontface(ZamEQ2_UI* ui) {
/*
	cairo_surface_t *bg;
	unsigned char * img_tmp;
	//img2surf((struct MyGimpImage const *) &img_eq2, &bg, &img_tmp);
	cairo_t *cr;
	if (!ui->frontface) {
		ui->frontface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, DAWIDTH, DAHEIGHT);
	}
	cr = cairo_create(ui->frontface);

	cairo_save(cr);
	cairo_set_source_surface(cr, bg, 0, 0);
	cairo_rectangle (cr, 0, 0, DAWIDTH, DAHEIGHT);
	cairo_fill(cr);
	cairo_paint(cr);
	cairo_restore(cr);

	cairo_destroy(cr);
	free(img_tmp);
*/

	cairo_t *cr;
	robtk_xydraw_set_surface(ui->xyp, NULL);
	ui->eqcurve = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, DAWIDTH, DAWIDTH);
	cr = cairo_create (ui->eqcurve);

	CairoSetSouerceRGBA(c_blk);
	cairo_set_operator (cr, CAIRO_OPERATOR_SOURCE);
	cairo_rectangle (cr, 0, 0, DAWIDTH, DAHEIGHT);
	cairo_fill (cr);

	cairo_save(cr);
	rounded_rectangle (cr, 10, 10, DAWIDTH - 20, DAHEIGHT - 20, 10);
	CairoSetSouerceRGBA(c_blk);
	cairo_fill_preserve(cr);
	cairo_clip(cr);
	cairo_set_operator (cr, CAIRO_OPERATOR_OVER);
	
	robtk_xydraw_set_surface(ui->xyp, ui->eqcurve);
}

static void calceqcurve(float val[], float x[], float y[])
{
	for (uint32_t i = 0; i < EQPOINTS; ++i) {
		x[i] = 0.0;
		y[i] = 0.0;
	}
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

static void xy_clip_fn(cairo_t *cr, void *data)
{
	rounded_rectangle(cr, 10, 10, DAWIDTH - 20, DAHEIGHT - 20, 10);
	cairo_clip(cr);
}


static bool cb_disp_changed (RobWidget* handle, void *data) {
	ZamEQ2_UI* ui = (ZamEQ2_UI*) (data);
	for (uint32_t i = 0; i < 4; ++i) {
		robwidget_show(ui->knob_freq[i]->rw, true);
		robwidget_show(ui->lbl_freq[i]->rw, true);
		robwidget_show(ui->knob_bw[i]->rw, true);
		robwidget_show(ui->lbl_bw[i]->rw, true);
		robwidget_show(ui->knob_gain[i]->rw, true);
		robwidget_show(ui->lbl_gain[i]->rw, true);
	}
	//robwidget_show(ui->knob_ingain->rw, true);
	robwidget_show(ui->knob_outgain->rw, true);
	robwidget_show(ui->lbl_ingain->rw, true);
	robwidget_show(ui->lbl_outgain->rw, true);
	
	return TRUE;
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

	return TRUE;
}

static void ui_disable(LV2UI_Handle handle)
{
}

static void ui_enable(LV2UI_Handle handle)
{
}

static void
size_request(RobWidget* handle, int *w, int *h) {
	*w = DAWIDTH;
	*h = DAHEIGHT;
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

	ui->xyp = robtk_xydraw_new(DAWIDTH, DAHEIGHT);
	//ui->xyp->rw->position_set = plot_position_right;
	robtk_xydraw_set_alignment(ui->xyp, 0, 0);
	robtk_xydraw_set_linewidth(ui->xyp, 1.5);
	robtk_xydraw_set_drawing_mode(ui->xyp, RobTkXY_ymax_zline);
	robtk_xydraw_set_mapping(ui->xyp, 1./EQPOINTS, 0, 1./EQPOINTS, 1.);
	robtk_xydraw_set_area(ui->xyp, 10, 10, DAWIDTH - 10, DAHEIGHT -10);
	robtk_xydraw_set_clip_callback(ui->xyp, xy_clip_fn, ui);
	robtk_xydraw_set_color(ui->xyp, 1.0, .0, .2, 1.0);

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
	
	robtk_xydraw_set_points(ui->xyp, EQPOINTS, ui->eqx, ui->eqy);

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

