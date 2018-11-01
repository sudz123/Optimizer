import './styles/main.scss';
import marked from 'marked';
import overview from './templates/overview';

const mainContent = document.getElementById('content');
mainContent.innerHTML = overview;
