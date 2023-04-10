#!/usr/bin/env k8

/*******************************
 * Command line option parsing *
 *******************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/************
 * Commands *
 ************/

function srf_cmd_enlong(args) {
	var c, opt = { min_len:200, dbl:false };
	while ((c = getopt(args, "l:d")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'd') opt.dbl = true;
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js enlong [options] <circ.fa>");
		print("Options:");
		print("  -l INT    min length [" + opt.min_len + "]");
		print("  -d        double sequences");
		return 1;
	}
	var buf = new Bytes();
	var file = new File(args[getopt.ind]);

	function enlong_seq(seq) {
		if (opt.dbl) return seq + "\n" + seq;
		var t = seq;
		while (t.length < opt.min_len)
			t += seq;
		return t;
	}

	var name = null, comm = null, seq = null;
	while (file.readline(buf) >= 0) {
		var m, l = buf.toString();
		if ((m = /^>(\S+)(.*)/.exec(l)) != null) {
			if (name != null) {
				print(">" + name + comm);
				print(enlong_seq(seq));
			}
			name = m[1], comm = m[2], seq = "";
		} else {
			seq += l;
		}
	}
	if (name != null) {
		print(">" + name + comm);
		print(enlong_seq(seq));
	}
	file.close();
	buf.destroy();
}

function srf_drop_paf(opt, t) {
	for (var i = 1; i <= 3; ++i) t[i] = parseInt(t[i]);
	for (var i = 6; i < t.length; ++i) t[i] = parseInt(t[i]);
	if (t[8] - t[7] < opt.min_len) return true;
	var clip = [t[7], t[6] - t[8]];
	if (t[4] == "+") {
		clip[0] = clip[0] < t[2]? clip[0] : t[2];
		clip[1] = clip[1] < t[1] - t[3]? clip[1] : t[1] - t[3];
	} else {
		clip[1] = clip[1] < t[2]? clip[1] : t[2];
		clip[0] = clip[0] < t[1] - t[3]? clip[0] : t[1] - t[3];
	}
	if (t[8] - t[7] < (t[8] - t[7] + clip[0] + clip[1]) * opt.min_cov) return true;
	var iden = t[9] / (t[10] + clip[0] + clip[1]);
	if (iden < opt.min_iden) return true;
	t[11] = (1 - iden) * 100;
	return false;
}

function srf_cmd_filter(args) {
	var c, opt = { min_len:190, min_cov:0.9, min_iden:0.0 };
	while ((c = getopt(args, "l:c:d:")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'c') opt.min_cov = parseFloat(getopt.arg);
		else if (c == 'd') opt.min_iden = parseFloat(getopt.arg);
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js filter [options] <in.paf>");
		print("Options:");
		print("  -l INT      min alignment length [" + opt.min_len + "]");
		print("  -c FLOAT    min coverage [" + opt.min_cov + "]");
		print("  -d FLOAT    min identity [" + opt.min_iden + "]");
		print("Notes: suggested minimap2 setting:");
		print("  minimap2 -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong circ.fa) reads.fa");
		return 1;
	}
	var buf = new Bytes();
	var file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t", 12);
		if (srf_drop_paf(opt, t)) continue;
		print(line);
	}
	file.close();
	buf.destroy();
}

function srf_cmd_bed2cnt(args) {
	var c, opt = { min_len:10000, thres2nd:0.5 };
	while ((c = getopt(args, "")) != null) {
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js bed2cnt <in.bed>");
		return 1;
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var hit = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		hit.push([t[3], t[0], parseInt(t[2]) - parseInt(t[1]), parseFloat(t[4])]);
	}
	file.close();
	buf.destroy();

	var a = hit.sort(function(x, y) { return x[0]<y[0]?-1:x[0]>y[0]?1:x[1]<y[1]?-1:x[1]>y[1]?1:0 });
	for (var i = 1, i0 = 0; i <= a.length; ++i) {
		if (i == a.length || a[i][0] != a[i0][0]) {
			var b = [];
			for (var j = i0 + 1, j0 = i0; j <= i; ++j) {
				if (j == i || a[j][1] != a[j0][1]) {
					var len = 0, score = 0;
					for (var k = j0; k < j; ++k)
						len += a[k][2], score += a[k][2] * a[k][3];
					var avg = score / len;
					score = len * Math.pow(1 - avg / 100, 4);
					if (len >= opt.min_len)
						b.push([a[j0][1], len, avg, score]);
					j0 = j;
				}
			}
			if (b.length > 0) {
				b = b.sort(function(x, y) { return y[3] - x[3]; });
				var j, thres = b[0][3] * opt.thres2nd;
				for (j = 0; j < b.length; ++j)
					if (b[j][3] < thres) break;
				b.length = j;
				var o = [];
				for (j = 1; j < b.length; ++j)
					o.push(b[j][0] + ":" + b[j][3].toFixed(2));
				print(a[i0][0], b[0][0], b[0][1], b[0][2].toFixed(2), b[0][3].toFixed(2), o.length, o.join("\t"));
			}
			i0 = i;
		}
	}
}

function srf_cmd_bed2abun(args) {
	var c, opt = { accumu:false, binplot:false, gsize:0, label:"", incl_all:false };
	while ((c = getopt(args, "abcg:l:")) != null) {
		if (c == 'a') opt.accumu = true;
		else if (c == 'b') opt.binplot = true;
		else if (c == 'g') opt.gsize = parseFloat(getopt.arg);
		else if (c == 'l') opt.label = getopt.arg;
		else if (c == 'c') opt.incl_all = true;
	}
	if (getopt.ind == args.length) {
		print("Usage: srfutils.js bed2abun [options] <in.bed>");
		print("Options:");
		print("  -g INT       total length [0]");
		print("  -a           generate data for accumulative plot");
		print("  -b           for stacked bar plot (<100,<1000,<2000,<5000,<10000,>=10000)");
		print("  -c           do not filter BED records by last column");
		print("  -l STR       label to output");
		return;
	}

	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var abun = {}, abun_all = {};
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (abun[t[3]] == null) abun[t[3]] = [0, 0];
		if (abun_all[t[3]] == null) abun_all[t[3]] = [0, 0];
		var len = parseInt(t[2]) - parseInt(t[1]);
		if (!opt.incl_all && parseInt(t[7]) > 0) {
			abun[t[3]][0] += len;
			abun[t[3]][1] += len * parseFloat(t[4]);
		}
		abun_all[t[3]][0] += len;
		abun_all[t[3]][1] += len * parseFloat(t[4]);
	}
	file.close();
	buf.destroy();
	var a = [];
	for (var x in abun)
		a.push([x, abun[x][0], abun[x][1] / abun[x][0], abun_all[x][0]]);
	var tot = 0;
	for (var i = 0; i < a.length; ++i)
		tot += a[i][1];
	if (opt.gsize > tot) tot = opt.gsize;
	if (opt.accumu || opt.binplot) {
		var b4 = [0, 0, 0, 0, 0, 0];
		for (var i = 0; i < a.length; ++i) {
			var m;
			if ((m = /-(\d+)$/.exec(a[i][0])) == null)
				throw Error("failed to match format");
			var len = parseInt(m[1]);
			a[i].push(len);
			if (len < 100) b4[0] += a[i][1];
			else if (len < 1000) b4[1] += a[i][1];
			else if (len < 2000) b4[2] += a[i][1];
			else if (len < 5000) b4[3] += a[i][1];
			else if (len < 10000) b4[4] += a[i][1];
			else b4[5] += a[i][1];
		}
		if (opt.accumu) {
			a.sort(function(x,y) { return x[3] - y[3] });
			for (var i = 0, c = 0; i < a.length; ++i) {
				c += a[i][1] / tot;
				print(c, a[i][3], a[i][1]);
			}
		} else {
			var label = opt.label != ""? opt.label : args[getopt.ind];
			for (var i = 0; i < b4.length; ++i)
				b4[i] /= tot;
			print(label, b4.join("\t"));
		}
	} else {
		a = a.sort(function(x,y) { return y[1] - x[1] });
		for (var i = 0; i < a.length; ++i)
			print(a[i][0], a[i][1], a[i][2], a[i][1] / tot, a[i][3] / tot);
	}
}

function srf_cmd_paf2bed(args) {
	var c, opt = { min_len:190, min_cov:0.9, min_iden:0.9, algo:1 };
	while ((c = getopt(args, "l:c:d:pa:")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'c') opt.min_cov = parseFloat(getopt.arg);
		else if (c == 'd') opt.min_iden = parseFloat(getopt.arg);
		else if (c == 'a') opt.algo = parseInt(getopt.arg);
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js paf2bed [options] <in.paf>");
		print("Options:");
		print("  -l INT      min alignment length [" + opt.min_len + "]");
		print("  -c FLOAT    min coverage [" + opt.min_cov + "]");
		print("  -d FLOAT    min identity [" + opt.min_iden + "]");
		print("  -a INT      algorithm (0: old; 1: new) [" + opt.algo + "]");
		print("Notes: suggested minimap2 setting:");
		print("  minimap2 -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong circ.fa) reads.fa");
		return 1;
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var b = [];

	function process0(b) {
		if (b.length == 0) return;
		b = b.sort(function(x,y) { return x[2]-y[2] });
		var tlen = b[0][1], j0 = 0;
		for (var j = 0; j < b.length; ++j) {
			var t = b[j];
			while (j0 < b.length && b[j0][3] <= t[2])
				++j0;
			var st = t[2], en = t[3];
			for (var i = j0; i < j; ++i) {
				var s = b[i];
				if (s[13] <= t[2]) continue;
				if (s[11] < t[11])
					st = st > s[13]? st : s[13];
			}
			for (var i = j0; i < j; ++i) {
				var s = b[i];
				if (s[13] <= st) continue;
				if (s[11] >= t[11])
					s[13] = st;
			}
			t[12] = st, t[13] = en;
		}
		var k = 0;
		for (var j = 0; j < b.length; ++j) {
			var t = b[j];
			if (t[13] - t[12] <= 0) continue;
			b[k++] = t;
		}
		b.length = k;
		if (b.length == 0) return;
		var en = b[0][13];
		var sc = (b[0][13] - b[0][12]) * b[0][11];
		var c = [];
		for (var j = 1, j0 = 0; j <= b.length; ++j) {
			if (j == b.length || b[j][0] != b[j0][0] || b[j][5] != b[j0][5] || b[j][12] != en) {
				//print(b[j0][0], b[j0][12], en, b[j0][5], sc / (en - b[j0][12]));
				var m, keep = 0, st = b[j0][12];
				if ((m = /-(\d+)$/.exec(b[j0][5])) == null)
					throw Error("Wrong contig format: " + b[j0][5]);
				var len = parseInt(m[1]);
				if (en - st > len * 2) keep = 1;
				else if ((st < 5 || tlen - en < 5) && en - st > len * 1.5) keep = 1;
				c.push([b[j0][0], st, en, b[j0][5], sc / (en - b[j0][12]), tlen, len, keep]);
				if (j < b.length)
					j0 = j, en = b[j][13], sc = (b[j][13] - b[j][12]) * b[j][11];
			} else {
				en = b[j][13], sc += (b[j][13] - b[j][12]) * b[j][11];
			}
		}
		b.length = 0;
		for (var j = 0; j < c.length; ++j) {
			if (!c[j][7]) {
				for (var k = j - 1; k >= 0; --k) {
					if (c[j][1] - c[k][2] >= c[j][6] * 1.5) break;
					if (c[j][3] == c[k][3]) c[j][7] = c[k][7] = 2;
				}
				for (var k = j + 1; k < c.length; ++k) {
					if (c[k][1] - c[j][2] >= c[j][6] * 1.5) break;
					if (c[j][3] == c[k][3]) c[j][7] = c[k][7] = 2;
				}
			}
			print(c[j].join("\t"));
		}
	}

	function output1(uj, v) {

		function appendv(v, uj, st, en) {
			var added = false;
			if (en == st) return;
			if (v.length > 0) {
				var ve = v[v.length - 1];
				if (ve[0] == uj[0] && ve[3] == uj[5] && ve[2] >= st && en > ve[1]) {
					if (ve[2] < en) { // something to add
						ve[4] = (ve[4] * (ve[2] - ve[1]) + uj[11] * (en - ve[2])) / (en - ve[1]); // update the weighted identity
						ve[2] = en;
					}
					added = true;
				}
			}
			if (!added) {
				var m, srf_len;
				srf_len = (m = /-(\d+)$/.exec(uj[5])) != null? parseInt(m[1]) : uj[6];
				v.push([uj[0], st, en, uj[5], uj[11], uj[1], srf_len, 0]);
			}
		}

		if (uj[14] == false) return;
		var w = uj[13], st1 = uj[12], x = [];
		for (var i = 0; i < w.length; ++i)
			if (w[i][1] > st1)
				x.push([w[i][0] > st1? w[i][0] : st1, w[i][1]]);
		if (x.length > 0) {
			var st = x[0][0], en = x[0][1], y = [];
			for (var i = 1; i < x.length; ++i) { // de-overlap blocked intervals
				if (x[i][0] > en) {
					y.push([st, en]);
					st = x[i][0], en = x[i][1];
				} else en = en > x[i][1]? en : x[i][1];
			}
			y.push([st, en]);
			// output unblocked intervals
			if (st1 < y[0][0])
				appendv(v, uj, st1, y[0][0]);
			for (var i = 1; i < y.length; ++i)
				appendv(v, uj, y[i-1][1], y[i][0]);
			if (y[y.length-1][1] < uj[3])
				appendv(v, uj, y[y.length-1][1], uj[3]);
		} else appendv(v, uj, st1, uj[3]);
		uj[14] = false;
	}

	function process1(b) {
		if (b.length == 0) return;
		b = b.sort(function(x,y) { return x[2]-y[2] });
		var u = [], v = [];
		while (b.length) {
			var t = b.shift(), st = t[2];
			for (var j = 0; j < u.length; ++j) {
				var uj = u[j];
				if (uj[0] != t[0] || uj[3] <= t[2]) { // t doesn't overlap with u[j]
					if (uj[14]) output1(uj, v);
				} else if (t[3] <= uj[3]) { // t is contained in u[j]
					if (t[11] >= uj[11]) st = t[3]; // u[j] is better
					else uj[13].push([t[2], t[3]]);
				} else { // t overlaps with u[j], but not contained
					if (t[11] >= uj[11]) st = st > uj[3]? st : uj[3];
					else uj[13].push([t[2], uj[3]]);
				}
			}
			var n_drop = 0;
			for (var j = 0; j < u.length; ++j) { // count out-of-range records from the beginning of u[]
				if (u[j][0] == t[0] && u[j][3] > t[2])
					break;
				else ++n_drop;
			}
			for (var j = 0; j < n_drop; ++j) // drop out-of-range records at the beginning of u[]
				u.shift();
			if (st >= t[3]) continue; // skip t
			t[12] = st, t[13] = [], t[14] = true; // t[13]: list of blocked intervals; t[14]: outputted or not
			u.push(t);
		}
		for (var j = 0; j < u.length; ++j)
			output1(u[j], v);
		v = v.sort(function(x,y) { return x[1]-y[1] });
		// filter by nearby hits
		for (var j = 0; j < v.length; ++j) {
			if (v[j][2] - v[j][1] >= v[j][6] * 2.0) v[j][7] = 1;
			else if (v[j][2] - v[j][1] >= v[j][6] * 1.5 && (v[j][5] - v[j][2] < 5 || v[j][1] < 5)) v[j][7] = 1; // relaxed at the edge
		}
		for (var j = 0; j < v.length; ++j) {
			var srf_len = v[j][6], lj = v[j][2] - v[j][1];
			if (v[j][7] == 1) continue;
			for (var k = j - 1; k >= 0; --k) {
				if (v[j][1] - v[k][2] > srf_len * 3.0) break;
				if (v[k][3] == v[j][3]) lj += v[k][2] - v[k][1];
			}
			for (var k = j + 1; k < v.length; ++k) {
				if (v[k][1] - v[j][2] > srf_len * 3.0) break;
				if (v[k][3] == v[j][3]) lj += v[k][2] - v[k][1];
			}
			if (lj >= srf_len * 2.0) {
				for (var k = j - 1; k >= 0; --k) {
					if (v[j][1] - v[k][2] > srf_len * 2.0) break;
					if (v[k][3] == v[j][3] && v[k][7] == 0) v[k][7] = 2;
				}
				for (var k = j + 1; k < v.length; ++k) {
					if (v[k][1] - v[j][2] > srf_len * 2.0) break;
					if (v[k][3] == v[j][3] && v[k][7] == 0) v[k][7] = 2;
				}
			}
		}
		// output
		for (var j = 0; j < v.length; ++j)
			print(v[j].join("\t"));
	}

	function process(b) {
		if (opt.algo == 0) process0(b);
		else if (opt.algo == 1) process1(b);
	}

	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t", 12);
		if (srf_drop_paf(opt, t)) continue;
		t[12] = t[2], t[13] = t[3];
		if (b.length > 0 && b[0][0] != t[0])
			process(b);
		b.push(t);
	}
	process(b);
	file.close();
	buf.destroy();
}

String.prototype.lpad = function(padString, length) {
	var str = this;
	while (str.length < length)
		str = padString + str;
	return str;
}

function srf_cmd_merge(args) {
	var c, opt = {};
	while ((c = getopt(args, "")) != null) {
	}
	if (args.length - getopt.ind < 2) {
		print("Usage: srfutils.js merge <prefix> <in.fa>");
		return;
	}
	var pre = args[getopt.ind];
	var buf = new Bytes();
	var file = new File(args[getopt.ind + 1]);
	var h = {}, n = 0, seq = '', name = null;
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if ((m = /^>(\S+)/.exec(line)) != null) {
			if (name != null) {
				if (h[seq] == null) {
					++n;
					print('>' + pre + '*' + String(n).lpad('0', 4) + ' ' + name);
					print(seq);
					h[seq] = 1;
				}
			}
			seq = '', name = m[1];
		} else seq += line;
	}
	if (name != null && h[seq] == null) {
		++n;
		print('>' + pre + '*' + String(n).lpad('0', 4) + ' ' + name);
		print(seq);
	}
	file.close();
	buf.destroy();
}

/*************************
 ***** main function *****
 *************************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: srfutils.js <command> [arguments]");
		print("Commands:");
		print("  enlong    enlong circ contigs");
		print("  filter    filter read-to-circ alignment");
		print("  paf2bed   extract non-overlapping regions with satellites");
		print("  bed2cnt   per-chromosome count");
		print("  bed2abun  per-contig abundance");
		print("  merge     merge identical sequences");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'filter') srf_cmd_filter(args);
	else if (cmd == 'enlong') srf_cmd_enlong(args);
	else if (cmd == 'merge') srf_cmd_merge(args);
	else if (cmd == 'paf2bed') srf_cmd_paf2bed(args);
	else if (cmd == 'bed2cnt') srf_cmd_bed2cnt(args);
	else if (cmd == 'bed2abun') srf_cmd_bed2abun(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
