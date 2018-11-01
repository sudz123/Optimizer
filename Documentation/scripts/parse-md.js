const path = require('path');
const fs = require('fs');
const docsDir = path.join(__dirname, '../docs');
const templatesDir = path.join(__dirname, '../src/templates');
const marked = require('marked');
const myMarked = marked;

myMarked.setOptions({
  renderer: new myMarked.Renderer(),
  highlight: code => require('highlight.js').highlightAuto(code).value,
  pedantic: false,
  gfm: true,
  tables: true,
  breaks: false,
  sanitize: false,
  smartLists: true,
  smartypants: false,
  xhtml: false,
});

for (const docPage of fs.readdirSync(docsDir)) {
  const docPath = path.join(docsDir, docPage);
  const doc = fs.readFileSync(docPath, 'utf-8');
  const parsed = myMarked(doc);
  const templateName = docPage.replace(/\..+/, '.js');
  const template = 'module.exports = `\n' + parsed + '\n`';
  // make template
  fs.writeFile(path.join(templatesDir, templateName), template, err => {
    if (err) {
      throw err;
    }
    console.log(`${templateName} created!`);
  });
}
