#!/usr/bin/env python
# -*- coding: utf-8 -*-


#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8')

from sys import argv
from jinja2 import Template
import re

example = """
# 宏基因组报告

## 项目整体流程概况

### 1.1 项目信息

&[](./images/general.txt)

### 1.2 分析项目及基本要求

| 分析项目                             | 分析要求      |
| :----------------------------------- | ------------- |
| 2.1 生物信息学分析                   |               |
| 2.1.1 测序原始数据的预处理和质控     | —             |
| 2.1.2 高质量序列的筛查过滤           | —             |
| 2.1.3 基于高质量序列的物种注释       | —             |
| 2.1.4 序列的组装拼接                 | —             |
| 2.1.5 非冗余序列集的构建             | —             |
| 2.1.6 基因预测                       | —             |
| 2.1.7 蛋白功能注释                   | —             |
| 2.2 功能组成分析                     |               |
| 2.2.1 各等级功能注释相对丰度分布分析 | —             |
| 2.2.2 相对丰度差异分析               | 样本（组）≥ 2 |
| 2.2.3 共有/独有KEGG代谢通路分析      | 样本（组）≥ 2 |
| 2.2.4 KEGG代谢通路富集分析           | 样本（组）≥ 2 |
| 2.2.5 NOG/CAZy功能比较分析           | 样本（组）≥ 2 |
| 2.2.6 功能注释丰度聚类分析           | 样本（组）≥ 2 |

...

## 项目分析结果

### 2.1 生物信息学分析

#### 2.1.1 测序原始数据的预处理和质控

本项目基于**Illumina NovaSeq**高通量测序平台...

![测序读长、插入片段和接头（Adapter）的构造图](./images/icon/adapter.png)

...

#### 2.1.3 基于高质量序列的物种注释

##### 2.1.3.1 kraken2

...
"""

html = """
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="renderer" content="webkit">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <title>{{ title }}</title>
        <link type="text/css" rel="stylesheet" href="static/css/bootstrap.min.css">
        <link type="text/css" rel="stylesheet" href="static/css/report.css">
        <link type="text/css" rel="stylesheet" href="static/css/bootstrap-table.css">
        <link type="text/css" rel="stylesheet" href="static/css/table.css">
        <script src="static/js/jquery-3.1.0.min.js" type="text/javascript"></script>
        <script src="static/js/bootstrap.min.js" type="text/javascript"></script>
        <script src="static/js/tables.js" type="text/javascript"></script>
        <script src="static/js/base.js" type="text/javascript"></script>
        <script src="static/js/bootstrap-table.js"></script>
        <script src="static/js/tableExport.js"></script>
        <script src="static/js/bootstrap-table-export.js"></script>
        <script src="static/js/tableList.js"></script>
    </head>
    <body data-spy="scroll" data-target=".scrollspy">
        <div id="slide-nav" class="d-print-none">
            <a href="http://www.personalbio.cn/" class="slide-nav-logo">
                <img src="./static/icon/logo.png" alt="" align="left">
            </a>
            <h2 class="slide-nav-title">{{ title }}</h2>
        </div>
        <div id="sidebar" class="fixed-column fixed-column-top bg-blue d-print-none bg-color">
            <ul class="nav justify-content-center">
            {% for page in pages %}
                {% if loop.first %}
                <li class="mb-3 justify-item{{ page.id }}">
                    <a class="active" data-toggle="pill" href="#part{{ page.id }}_tab_pane">{{ page.title }}</a>
                </li>
                {% else %}
                <li class="mb-3 justify-item{{ page.id }}">
                    <a class="" data-toggle="pill" href="#part{{ page.id }}_tab_pane">{{ page.title }}</a>
                </li>
                {% endif %}
            {% endfor %}
            </ul>
        </div>
        <div id="frontpage" class="d-none d-print-block">
            <img src="static/icon/page_bg.jpg">
            <div id="pageTitle">
                <p>项目名称：{{ title }}</p>
                <p>委托单位：XXXXXXXX</p>
                <p>制定日期：YYYYYYYY</p>
            </div>
        </div>
        
        <div id="floatbar" class="d-print-none">
            <div id="floatbar_goTop">
                <span class="glyphicon glyphicon-chevron-up"></span>
                <p>返回顶部</p>
            </div>
            <div class="floatbar-wrap">
                <div class="floatbar-file float-mouse">
                    <a href="">
                        <span class="glyphicon glyphicon-save-file"></span>
                        <p>文件下载</p>
                    </a>
                </div>
                <div class="float_btn float-mouse">
                    <span class="glyphicon glyphicon-envelope"></span>
                    <p>联系方式</p>
                </div>
                <!-- <div class="float_btn"><span class="glyphicon glyphicon-qrcode"></span><p>微信订阅号</p></div> -->
                <div class="floatbar-info float-mouse">
                    <a id="fankui" href="http://c.eqxiu.com/s/U1PFE5Pk" target="_blank">
                        <span class="glyphicon glyphicon-info-sign"></span>
                        <p>信息反馈</p>
                    </a>
                </div>
                <!-- <div id="floatbar_remove"><span class="glyphicon glyphicon-remove"></span></div> -->
                <!-- <div id="floatbar_goTop"><span class="glyphicon glyphicon-chevron-up"></span><p>返回顶部</p></div> -->
            </div>
            <div class="packup packup-active">
                <i class="glyphicon glyphicon-plus"></i>
            </div>
        </div>
        <div class="float_tips tel-box">
            <div class="tel-left">
                <p class="tel-left-p">
                    <b>版权所有：</b>
                    <br>派森诺生物单细胞空转产品部
                </p>
                <p class="tel-left-p">
                    <b>邮箱：</b>
                    <br>sc_support@personalbio.cn
                </p>
            </div>
            <div class="tel-left-border">
                <div></div>
            </div>
            <img src="static/icon/2Dplot_1.jpg" alt="" class="tel-right">
        </div>
        <div class="float_tips">
            <img src="static/icon/2Dplot_1.jpg" alt="">
        </div>
        <div class="float_tips">
            <img src="static/icon/2Dplot_2.png" alt="">
        </div>
        <div class="tab-content">
            {% for page in pages %}
            <div id="part{{ page.id }}_tab_pane" class="tab-pane {% if loop.first %}active {% endif %}bg-gray d-print-block">
                <div class="catalog fixed-column bg-white d-print-none catalog-id1">
                    <!-- <a href="http://www.personalbio.cn" target="_blank" class="tab-content-a"><img src="static/icon/logo.png" link="www.personalbio.cn"></a> -->
                    <nav class="nav-pills scrollspy nav-pills-top nav-pills-one">
                    {% for part in page.parts %}
                        <a class="nav-link {% if loop.first %}active-color{% endif %}" href="#part{{ page.id }}_ch{{ part.id }}">{{ page.id }}.{{ part.id }} &nbsp;{{ part.title }}</a>
                        {% if part.subparts|length > 0 %}
                        <nav class="nav-pills nav-pills-two">
                            {% for subpart in part.subparts %}
                            <a class="nav-link" href="#part{{ page.id }}_ch{{ part.id }}_{{ subpart.id }}">{{ page.id }}.{{ part.id }}.{{ subpart.id }} &nbsp;{{ subpart.title }}</a>
                            {% endfor %}
                        </nav>
                        {% endif %}
                    {% endfor %}
                    </nav>
                </div>
                <div class="content bg-white content-id{{ page.id }}">
                    <div class="col-12 border-bottom mb-5 d-none d-print-block" style="page-break-before: always;">
                        <img src="static/icon/logo.png" class="col-2 offset-10 mb-2"/>
                    </div>
                    <div class="partSeg row text-center d-none d-print-block" style="margin-top: 40%">
                        <h1 class="display-4 font-weight-bold">{{ page.title }}</h1>
                    </div>
                    {% set outer_loop = loop %}
                    {% for part in page.parts %}
                    <div id="part{{ page.id }}_ch{{ part.id }}" class="mar-bot-80">
                        {% if not outer_loop.first or loop.index != 2 %}
                        <div class="col-12 border-bottom mb-5 d-none d-print-block" style="page-break-before: always;">
                            <img src="static/icon/logo.png" class="col-2 offset-10 mb-2"/>
                        </div>
                        {% endif %}
                        <h2 class="chapter-title base-item1">{{ page.id }}.{{ part.id }} &nbsp;{{ part.title }}</h2>
                        {% if part.subparts|length > 0 %}
                        {% for subpart in part.subparts %}
                        <div id="part{{ page.id }}_ch{{ part.id }}_{{ subpart.id }}" class="mar-bot-80">
                            <h3 class="section-title base-item10">{{ page.id }}.{{ part.id }}.{{ subpart.id }} &nbsp;{{ subpart.title }}</h3>
                            {{ subpart.text }}
                        </div>
                        {% endfor %}
                        {% else %}
                        {{ part.text }}
                        {% endif %}
                    </div>
                    {% endfor %}
                </div>
            </div>
            {% endfor %}
        </div>
        <div class="carousel-img-zoom">
            <img src="./result/img/1.png" alt="">
        </div>
        <script src="static/js/main.js"></script>
    </body>
</html>
"""

para = """
<p class="paragraph">{text}</p>
"""

bold = """<strong>{text}</strong>"""

link = """<a href="{url}" target="_blank">{name}</a>"""

img = """
                            <div component_type="img">
                                <div class="row">
                                    <div class="col-xl-12 col-lg-12 col-12 offset-xl-0 offset-lg-0 offset-0">
                                        <img class="img-fluid" alt="" src="{url}">
                                        <b class="imgname col-12">{title}</b>
                                        <p class="tablenote col-12">
                                        {text}
                                        </p>
                                    </div>
                                </div>
                            </div>
"""

table = """
                            <div class="row" component_type="table">
                                <div class="col-xl-10 col-lg-10 col-10 offset-xl-1 offset-lg-1 offset-1">
                                    <table class="onehead_table three-line-table table-sm" id="">
                                        <tbody>
                                        {% for item in table %}
                                        <tr>
                                            {% for it in item %}
                                            <th>{{ it }}</th>
                                            {% endfor %}
                                        </tr>
                                        {% endfor %}
                                        </tbody>
                                    </table>
                                    </div>
                            </div>
"""

table_link = """
                            <div class="row" component_type="table">
                                <div class="col-xl-10 col-lg-10 col-12 offset-xl-1 offset-lg-1">
                                    <b class="tablename col-12">{{ title }}</b>
                                    <div id="reportTableDiv">
                                    <table id="reportTable{{ num }}"></table>
                                    </div>
                                    <script type="text/javascript">
                                        $('#reportTable{{ num }}').bootstrapTable({
                                            height: 300,
                                            pagination: true,
                                            pageSize: 10,
                                            pageNumber: 1,
                                            showColumns: true,
                                            showExport: true,
                                            exportDataType: 'basic',
                                            exportTypes: ['csv', 'xls', 'doc', 'txt'],
                                            search: true,
                                            tableTitle: '{{ table.title }}',
                                            optionalParams: {
                                                fileName: '{{ table.title }}',
                                            },
                                            titleContent: {{ table.content }},
                                            columns: [
                                            {% for col in table.cols %}
                                            {field:"{{ col }}",title:"{{ col }}",align:"center",valign:"middle",sortable:"true"},
                                            {% endfor %}
                                            ],
                                            data: [
                                            {% for data in table.datas %}
                                            {{ data.__str__() }},
                                            {% endfor %}
                                            ],
                                            ellipsis: false,
                                            formatLoadingMessage: function () {
                                                return "请稍等，正在加载中...";
                                            },
                                            onSearch: function (text) {
                                            },
                                            onPageChange: function (size, number) {

                                            },
                                            formatNoMatches: function() {
                                                return '暂无相关内容！';
                                            }
                                        });
                                        $(window).resize(function () {
                                            $('#reportTable{{ num }}').bootstrapTable('resetView');
                                        });
                                    </script>
                                </div>
                            </div>
"""

def parse_table(item):
    with open(item) as t:
        content = t.read().strip().split("\n")

    content = [ line.split("\t") for line in content ]
    data = [ dict(zip(content[0], line)) for line in content[1:] ]
    table = {
        "title" : "",
        "cols" : content[0],
        "content" : content[0].__str__(),
        "datas" : data
    }
    return table

def insert_table(tbl):
    global table
    t = Template(table)
    return t.render(table=tbl)

def cons_page(page):
    pages = []
    for p in page:
        if p["level"] == "h2":
            pages.append(p)
        elif p["level"] == "h3":
            pages[-1]["parts"].append(p)
        elif p["level"] == "h4":
            pages[-1]["parts"][-1]["subparts"].append(p)
            
    return pages
    
    
def parse(md_f):
    html_text = ""
    global html
    page = []
    table = []
    table_cnt = 0
    id = 0
    num = 0
    subnum = 0
    text = ""
    with open(md_f) as md:
        for line in md:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith("|"):
                if not "----" in line:
                    row = line.split("|")
                    row = [ item.strip() for item in row ]
                    table.append(row[1:-1])
            else:
                if table:
                    text += insert_table(table)
                    table = []
                
            row = line.split(" ")
            title = " ".join(row[1:])

            if row[0] == "#":
                main_title = title
            elif row[0] == "##":
                if page:
                    page[-1]["text"] = text
                    text = ""
                id += 1
                num = 0
                subnum = 0
                page.append({
                    "level" : "h2",
                    "id" : id,
                    "title" : title,
                    "parts" : []
                })
                
            elif row[0] == "###":
                if page:
                    page[-1]["text"] = text
                    text = ""
                num += 1
                subnum = 0
                page.append({
                    "level" : "h3",
                    "id" : num,
                    "title" : title,
                    "text" : "",
                    "subparts" : []
                })
            elif row[0] == "####":
                if page:
                    page[-1]["text"] = text
                    text = ""
                subnum += 1
                page.append({
                    "level" : "h4",
                    "id" : subnum,
                    "title" : title,
                    "text" : "",
                })
            elif row[0] == "#####":
                text += '<h4 class="chapter-title base-item8">{}</h4>'.format(title.replace(' ', ' &nbsp;'))
            elif line.startswith("!"):
                name = line[2:line.find("]")]
                url = line[line.find("(")+1:line.find(")")]
                global img
                text += img.format(url=url, title=name, text="")
            elif line.startswith("&"):
                name = line[2:line.find("]")]
                url = line[line.find("(")+1:line.find(")")]
                table_cnt += 1
                global table_link
                t = Template(table_link)
                text += t.render(table=parse_table(url), title=name, num=table_cnt)
            elif line.startswith("|"):
                pass
            elif line.startswith("<"):
                text += line
            elif line.endswith(">"):
                text += line
            else:
                # bold
                strong = re.compile('\*\*([^\*]*)\*\*')
                image = re.compile('\[([^\[]*)\]\(([^\(]*)\)')
                italic = re.compile('\*([^\*]*)\*')
                line = re.sub(strong, lambda x: bold.format(text=x.groups()[0]), line)
                line = re.sub(image, lambda x: link.format(name=x.groups()[0], url=x.groups()[1]), line)
                line = re.sub(italic, lambda x: '<i>{}</i>'.format(x.groups()[0]), line)
                text += para.format(text=line)
                
    page[-1]["text"] = text
    pages = cons_page(page)
    global html
    t = Template(html)
    return t.render(pages=pages, title=main_title)
    
def main():
    print(parse(argv[1]))
    

if __name__ == "__main__":
    main()
