---
layout: page
title: "Undergraduate Research"
permalink: /r-undergrad/
---

<h1>Undergraduate Research</h1>

<ul>
  {% for post in site.posts %}
    {% if post.categories contains "r-undergrad" %}
      <li>
        <a href="{{ post.url }}">{{ post.title }}</a> - {{ post.date | date: "%B %d, %Y" }}
      </li>
    {% endif %}
  {% endfor %}
</ul>
