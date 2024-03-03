import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';

import { AuthModule } from '@auth0/auth0-angular';
import { ButtonModule } from 'primeng/button';
import {MenubarModule} from 'primeng/menubar';
import {InputTextModule} from 'primeng/inputtext'
import {CardModule} from 'primeng/card';
import {ReactiveFormsModule} from '@angular/forms'; 
import {DropdownModule} from 'primeng/dropdown';

import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { HomeComponent } from './home/home.component';
import { NavbarComponent } from './navbar/navbar.component';
import { AuthButtonComponent } from './auth-button/auth-button.component';
import { UserFormComponent } from './user-form/user-form.component';

@NgModule({
  declarations: [
    AppComponent,
    HomeComponent,
    NavbarComponent,
    AuthButtonComponent,
    UserFormComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    ButtonModule,
    MenubarModule,
    CardModule,
    InputTextModule,
    DropdownModule,
    ReactiveFormsModule,
    AuthModule.forRoot({
      domain: 'dev-jz7lwt5lmlxvurxt.us.auth0.com',
      clientId: 'FPoqV6Tf7mvnhroEZ5wqSOhbrs3G0P2k',
      authorizationParams: {
        redirect_uri: window.location.origin
      }
    }),
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
